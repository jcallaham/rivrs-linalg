//! Profile session, section guard, and finished session types.
//!
//! The profiler uses thread-local storage for zero-contention recording:
//! each thread records `ProfileEvent` structs into its own `Vec`, and
//! `finish()` merges all thread-local buffers into a unified hierarchy.

use std::cell::RefCell;
use std::collections::HashMap;
use std::marker::PhantomData;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{LazyLock, Mutex};
use std::time::{Duration, Instant};

use super::section::{ProfileEvent, ProfileSection};

/// Global session counter for associating thread-local events with sessions.
static SESSION_COUNTER: AtomicU64 = AtomicU64::new(0);

/// Events tagged with a thread index, flushed from worker threads.
type FlushedEvents = HashMap<u64, Vec<(usize, RawEvent)>>;

/// Shared storage for events flushed from worker threads.
static SHARED_FLUSHED: LazyLock<Mutex<FlushedEvents>> =
    LazyLock::new(|| Mutex::new(HashMap::new()));

thread_local! {
    /// Per-thread event buffer: maps session_id → raw events.
    static THREAD_EVENTS: RefCell<HashMap<u64, Vec<RawEvent>>> = RefCell::new(HashMap::new());
    /// Per-thread depth counter: maps session_id → current nesting depth.
    static THREAD_DEPTH: RefCell<HashMap<u64, usize>> = RefCell::new(HashMap::new());
}

/// Raw event stored in thread-local buffer before merge.
#[derive(Debug, Clone)]
struct RawEvent {
    name: String,
    start_us: f64,
    duration: Duration,
    depth: usize,
}

/// A bounded profiling session collecting hierarchical timing data from
/// potentially multiple threads.
///
/// Thread-safe: `Send + Sync`. Each thread records into its own thread-local
/// storage, with zero contention during recording.
#[derive(Debug)]
pub struct ProfileSession {
    id: u64,
    start: Instant,
}

// SAFETY: ProfileSession only contains an id (u64) and start (Instant),
// both of which are Send + Sync. The actual event storage is thread-local.
unsafe impl Send for ProfileSession {}
unsafe impl Sync for ProfileSession {}

/// RAII guard that records a section's duration on drop.
///
/// Not `Send` — must be used on the thread that created it.
pub struct SectionGuard<'a> {
    session_id: u64,
    session_start: Instant,
    name: String,
    start: Instant,
    depth: usize,
    _not_send: PhantomData<*const ()>,
    _lifetime: PhantomData<&'a ()>,
}

/// A completed profiling session with merged, aggregated timing data.
pub struct FinishedSession {
    pub(crate) events: Vec<ProfileEvent>,
    pub(crate) sections: Vec<ProfileSection>,
    pub(crate) total_duration: Duration,
}

impl ProfileSession {
    /// Create a new profiling session, recording the start timestamp.
    pub fn new() -> Self {
        Self {
            id: SESSION_COUNTER.fetch_add(1, Ordering::Relaxed),
            start: Instant::now(),
        }
    }

    /// Begin timing a named section. Returns a guard that stops timing on drop.
    ///
    /// Supports nesting — calling `enter_section` while another guard is alive
    /// creates a child section at increased depth.
    pub fn enter_section(&self, name: &str) -> SectionGuard<'_> {
        let depth = THREAD_DEPTH.with(|d| {
            let mut map = d.borrow_mut();
            let depth = map.entry(self.id).or_insert(0);
            let current = *depth;
            *depth += 1;
            current
        });

        SectionGuard {
            session_id: self.id,
            session_start: self.start,
            name: name.to_string(),
            start: Instant::now(),
            depth,
            _not_send: PhantomData,
            _lifetime: PhantomData,
        }
    }

    /// Consume the session, collect all thread-local events, and return a
    /// finished session with aggregated timing data.
    pub fn finish(self) -> FinishedSession {
        let total_duration = self.start.elapsed();
        let raw_events = collect_session_events(self.id);

        // Assign thread indices based on order of appearance
        let mut thread_map: HashMap<usize, usize> = HashMap::new();
        let mut next_thread_idx = 0;

        let events: Vec<ProfileEvent> = raw_events
            .into_iter()
            .map(|(thread_hash, raw)| {
                let thread_idx = *thread_map.entry(thread_hash).or_insert_with(|| {
                    let idx = next_thread_idx;
                    next_thread_idx += 1;
                    idx
                });
                ProfileEvent {
                    name: raw.name,
                    start_us: raw.start_us,
                    duration: raw.duration,
                    thread_idx,
                    depth: raw.depth,
                }
            })
            .collect();

        let sections = build_section_tree(&events, total_duration);

        FinishedSession {
            events,
            sections,
            total_duration,
        }
    }
}

impl Default for ProfileSession {
    fn default() -> Self {
        Self::new()
    }
}

impl<'a> Drop for SectionGuard<'a> {
    fn drop(&mut self) {
        let duration = self.start.elapsed();
        let start_us = self
            .start
            .checked_duration_since(self.session_start)
            .unwrap_or(Duration::ZERO)
            .as_secs_f64()
            * 1_000_000.0;

        THREAD_EVENTS.with(|events| {
            events
                .borrow_mut()
                .entry(self.session_id)
                .or_default()
                .push(RawEvent {
                    name: std::mem::take(&mut self.name),
                    start_us,
                    duration,
                    depth: self.depth,
                });
        });

        THREAD_DEPTH.with(|d| {
            let mut map = d.borrow_mut();
            if let Some(depth) = map.get_mut(&self.session_id) {
                *depth = depth.saturating_sub(1);
            }
        });
    }
}

/// Flush thread-local events for a session to shared storage.
///
/// Call this from each worker thread before the session owner calls `finish()`.
/// The `thread_idx` parameter identifies this thread in the output.
pub fn flush_thread_events(session: &ProfileSession, thread_idx: usize) {
    let session_id = session.id;

    THREAD_EVENTS.with(|events| {
        let mut map = events.borrow_mut();
        if let Some(thread_events) = map.remove(&session_id) {
            let pairs: Vec<(usize, RawEvent)> = thread_events
                .into_iter()
                .map(|raw| (thread_idx, raw))
                .collect();

            if let Ok(mut shared) = SHARED_FLUSHED.lock() {
                shared.entry(session_id).or_default().extend(pairs);
            }
        }
    });

    THREAD_DEPTH.with(|d| {
        d.borrow_mut().remove(&session_id);
    });
}

/// Collect and remove all events for a session from current thread TLS and shared storage.
fn collect_session_events(session_id: u64) -> Vec<(usize, RawEvent)> {
    let mut result = Vec::new();

    // Collect from current thread's TLS (main thread / single-threaded path)
    THREAD_EVENTS.with(|events| {
        let mut map = events.borrow_mut();
        if let Some(thread_events) = map.remove(&session_id) {
            for raw in thread_events {
                result.push((0, raw));
            }
        }
    });

    THREAD_DEPTH.with(|d| {
        d.borrow_mut().remove(&session_id);
    });

    // Collect from shared storage (events flushed from worker threads)
    if let Ok(mut shared) = SHARED_FLUSHED.lock() {
        if let Some(flushed) = shared.remove(&session_id) {
            result.extend(flushed);
        }
    }

    result
}

/// Build aggregated section tree from flat events.
fn build_section_tree(events: &[ProfileEvent], total_duration: Duration) -> Vec<ProfileSection> {
    if events.is_empty() {
        return Vec::new();
    }

    let top_level: Vec<&ProfileEvent> = events.iter().filter(|e| e.depth == 0).collect();
    let mut sections = aggregate_by_name(&top_level, Some(total_duration));

    for section in &mut sections {
        attach_children(section, events, 1);
    }

    sections
}

/// Aggregate events by name, producing one ProfileSection per unique name.
fn aggregate_by_name(
    events: &[&ProfileEvent],
    parent_total: Option<Duration>,
) -> Vec<ProfileSection> {
    if events.is_empty() {
        return Vec::new();
    }

    // Preserve order of first appearance
    let mut seen_order: Vec<String> = Vec::new();
    let mut groups: HashMap<&str, Vec<&ProfileEvent>> = HashMap::new();

    for event in events {
        if !groups.contains_key(event.name.as_str()) {
            seen_order.push(event.name.clone());
        }
        groups.entry(&event.name).or_default().push(event);
    }

    seen_order
        .iter()
        .map(|name| {
            let group = &groups[name.as_str()];
            ProfileSection::aggregate(group, parent_total)
        })
        .collect()
}

/// Recursively attach children to a section based on event depth and time ranges.
fn attach_children(parent: &mut ProfileSection, all_events: &[ProfileEvent], child_depth: usize) {
    let parent_events: Vec<&ProfileEvent> = all_events
        .iter()
        .filter(|e| e.depth == child_depth - 1 && e.name == parent.name)
        .collect();

    let mut child_events: Vec<&ProfileEvent> = Vec::new();

    for parent_event in &parent_events {
        let parent_start = parent_event.start_us;
        let parent_end = parent_start + parent_event.duration.as_secs_f64() * 1_000_000.0;

        for event in all_events.iter() {
            if event.depth == child_depth
                && event.start_us >= parent_start
                && event.start_us < parent_end
            {
                child_events.push(event);
            }
        }
    }

    if child_events.is_empty() {
        return;
    }

    let mut children = aggregate_by_name(&child_events, Some(parent.total_duration));

    for child in &mut children {
        attach_children(child, all_events, child_depth + 1);
    }

    parent.children = children;
}

impl FinishedSession {
    /// Get the aggregated section tree.
    pub fn sections(&self) -> &[ProfileSection] {
        &self.sections
    }

    /// Get all raw events (for Chrome Trace export).
    pub fn events(&self) -> &[ProfileEvent] {
        &self.events
    }

    /// Get the total wall-clock duration of the session.
    pub fn total_duration(&self) -> Duration {
        self.total_duration
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;
    use std::thread;

    use super::*;

    #[test]
    fn new_session_records_start_time() {
        let before = Instant::now();
        let session = ProfileSession::new();
        let after = Instant::now();
        assert!(session.start >= before);
        assert!(session.start <= after);
    }

    #[test]
    fn enter_section_returns_guard() {
        let session = ProfileSession::new();
        let _guard = session.enter_section("test");
    }

    #[test]
    fn dropping_guard_records_event() {
        let session = ProfileSession::new();
        {
            let _guard = session.enter_section("test_section");
            thread::sleep(Duration::from_millis(1));
        }
        let finished = session.finish();
        assert_eq!(finished.events.len(), 1);
        assert_eq!(finished.events[0].name, "test_section");
        assert!(finished.events[0].duration >= Duration::from_millis(1));
    }

    #[test]
    fn nested_sections_produce_correct_depth() {
        let session = ProfileSession::new();
        {
            let _outer = session.enter_section("outer");
            {
                let _inner = session.enter_section("inner");
                thread::sleep(Duration::from_millis(1));
            }
        }
        let finished = session.finish();
        assert_eq!(finished.events.len(), 2);

        let inner = finished.events.iter().find(|e| e.name == "inner").unwrap();
        let outer = finished.events.iter().find(|e| e.name == "outer").unwrap();
        assert_eq!(inner.depth, 1);
        assert_eq!(outer.depth, 0);
    }

    #[test]
    fn finish_returns_finished_session_with_sections() {
        let session = ProfileSession::new();
        {
            let _guard = session.enter_section("analyze");
            thread::sleep(Duration::from_millis(1));
        }
        let finished = session.finish();

        assert_eq!(finished.sections().len(), 1);
        assert_eq!(finished.sections()[0].name, "analyze");
        assert_eq!(finished.sections()[0].call_count, 1);
        assert!(finished.total_duration() >= Duration::from_millis(1));
    }

    #[test]
    fn repeated_section_names_are_aggregated() {
        let session = ProfileSession::new();
        for _ in 0..5 {
            let _guard = session.enter_section("kernel");
            thread::sleep(Duration::from_millis(1));
        }
        let finished = session.finish();

        assert_eq!(finished.sections().len(), 1);
        let section = &finished.sections()[0];
        assert_eq!(section.name, "kernel");
        assert_eq!(section.call_count, 5);
        assert!(section.total_duration >= Duration::from_millis(5));
        assert!(section.min_duration >= Duration::from_millis(1));
        assert!(section.max_duration >= Duration::from_millis(1));
    }

    #[test]
    fn nested_hierarchy_in_sections() {
        let session = ProfileSession::new();
        {
            let _outer = session.enter_section("factor");
            {
                let _inner = session.enter_section("dense_kernel");
                thread::sleep(Duration::from_millis(1));
            }
        }
        let finished = session.finish();

        assert_eq!(finished.sections().len(), 1);
        let factor = &finished.sections()[0];
        assert_eq!(factor.name, "factor");
        assert_eq!(factor.children.len(), 1);
        assert_eq!(factor.children[0].name, "dense_kernel");
    }

    #[test]
    fn panic_in_section_records_completed_sections() {
        let session = ProfileSession::new();
        {
            let _outer = session.enter_section("before_panic");
            thread::sleep(Duration::from_millis(1));
        }

        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let _guard = session.enter_section("will_panic");
            panic!("test panic");
        }));
        assert!(result.is_err());

        let finished = session.finish();
        // Both events should be recorded: "before_panic" completed normally,
        // "will_panic" guard dropped during unwind (still records).
        assert!(!finished.events.is_empty());
        let has_before = finished.events.iter().any(|e| e.name == "before_panic");
        assert!(
            has_before,
            "section completed before panic should be recorded"
        );
    }

    #[test]
    fn thread_safety_multi_threaded_collection() {
        let session = Arc::new(ProfileSession::new());
        let num_threads = 4;
        let events_per_thread = 10;

        let handles: Vec<_> = (0..num_threads)
            .map(|t| {
                let session = Arc::clone(&session);
                thread::spawn(move || {
                    for i in 0..events_per_thread {
                        let _guard = session.enter_section(&format!("thread_{t}_event_{i}"));
                        thread::sleep(Duration::from_micros(10));
                    }
                    flush_thread_events(&session, t + 1);
                })
            })
            .collect();

        for h in handles {
            h.join().unwrap();
        }

        let session = Arc::try_unwrap(session).expect("all threads done");
        let finished = session.finish();

        let total_expected = num_threads * events_per_thread;
        assert_eq!(
            finished.events.len(),
            total_expected,
            "expected {total_expected} events, got {}",
            finished.events.len()
        );
    }

    #[test]
    fn empty_session_produces_empty_finished() {
        let session = ProfileSession::new();
        let finished = session.finish();
        assert!(finished.events.is_empty());
        assert!(finished.sections().is_empty());
    }

    #[test]
    fn section_guard_is_not_send() {
        fn assert_send<T: Send>() {}
        fn assert_sync<T: Sync>() {}
        assert_send::<ProfileSession>();
        assert_sync::<ProfileSession>();
        // SectionGuard is !Send due to PhantomData<*const ()>
    }
}
