//! Profile event and section data types.

use std::time::Duration;

/// A raw timing record from a single section entry/exit on a single thread.
#[derive(Debug, Clone)]
pub struct ProfileEvent {
    /// Section label.
    pub name: String,
    /// Absolute start time (microseconds from session start).
    pub start_us: f64,
    /// Duration of this section invocation.
    pub duration: Duration,
    /// Thread index (0-based, assigned during merge).
    pub thread_idx: usize,
    /// Nesting depth (0 = top level).
    pub depth: usize,
}

/// An aggregated view of one or more invocations of a named section.
#[derive(Debug, Clone)]
pub struct ProfileSection {
    /// Section label.
    pub name: String,
    /// Sum of all invocation durations.
    pub total_duration: Duration,
    /// Number of times this section was entered.
    pub call_count: u64,
    /// Shortest invocation.
    pub min_duration: Duration,
    /// Longest invocation.
    pub max_duration: Duration,
    /// Average duration (total / count).
    pub mean_duration: Duration,
    /// Percentage of parent section's total time (0.0-100.0).
    pub parent_pct: f64,
    /// Nested child sections.
    pub children: Vec<ProfileSection>,
}

impl ProfileSection {
    /// Aggregate a group of events with the same name and depth into a section.
    ///
    /// `parent_total` is the parent section's total duration, used to compute `parent_pct`.
    /// If `None`, `parent_pct` is set to 0.0 (top-level sections use session total later).
    pub(crate) fn aggregate(events: &[&ProfileEvent], parent_total: Option<Duration>) -> Self {
        assert!(!events.is_empty(), "cannot aggregate empty event list");

        let name = events[0].name.clone();
        let total_duration: Duration = events.iter().map(|e| e.duration).sum();
        let call_count = events.len() as u64;
        let min_duration = events.iter().map(|e| e.duration).min().unwrap();
        let max_duration = events.iter().map(|e| e.duration).max().unwrap();
        let mean_duration = total_duration / call_count as u32;

        let parent_pct = match parent_total {
            Some(pt) if !pt.is_zero() => total_duration.as_secs_f64() / pt.as_secs_f64() * 100.0,
            _ => 0.0,
        };

        Self {
            name,
            total_duration,
            call_count,
            min_duration,
            max_duration,
            mean_duration,
            parent_pct,
            children: Vec::new(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_event(name: &str, start_us: f64, dur_ms: u64, depth: usize) -> ProfileEvent {
        ProfileEvent {
            name: name.to_string(),
            start_us,
            duration: Duration::from_millis(dur_ms),
            thread_idx: 0,
            depth,
        }
    }

    #[test]
    fn profile_event_construction_and_fields() {
        let event = make_event("analyze", 100.0, 50, 0);
        assert_eq!(event.name, "analyze");
        assert_eq!(event.start_us, 100.0);
        assert_eq!(event.duration, Duration::from_millis(50));
        assert_eq!(event.thread_idx, 0);
        assert_eq!(event.depth, 0);
    }

    #[test]
    fn profile_event_clone() {
        let event = make_event("factor", 0.0, 100, 1);
        let cloned = event.clone();
        assert_eq!(cloned.name, event.name);
        assert_eq!(cloned.duration, event.duration);
    }

    #[test]
    fn aggregate_single_event() {
        let e = make_event("analyze", 0.0, 100, 0);
        let section = ProfileSection::aggregate(&[&e], None);
        assert_eq!(section.name, "analyze");
        assert_eq!(section.total_duration, Duration::from_millis(100));
        assert_eq!(section.call_count, 1);
        assert_eq!(section.min_duration, Duration::from_millis(100));
        assert_eq!(section.max_duration, Duration::from_millis(100));
        assert_eq!(section.mean_duration, Duration::from_millis(100));
        assert!(section.children.is_empty());
    }

    #[test]
    fn aggregate_repeated_labels() {
        let events = [
            make_event("kernel", 0.0, 10, 1),
            make_event("kernel", 100.0, 20, 1),
            make_event("kernel", 200.0, 30, 1),
        ];
        let refs: Vec<&ProfileEvent> = events.iter().collect();
        let parent_total = Duration::from_millis(100);
        let section = ProfileSection::aggregate(&refs, Some(parent_total));

        assert_eq!(section.name, "kernel");
        assert_eq!(section.call_count, 3);
        assert_eq!(section.total_duration, Duration::from_millis(60));
        assert_eq!(section.min_duration, Duration::from_millis(10));
        assert_eq!(section.max_duration, Duration::from_millis(30));
        assert_eq!(section.mean_duration, Duration::from_millis(20));
        // 60ms / 100ms = 60%
        assert!((section.parent_pct - 60.0).abs() < 0.1);
    }

    #[test]
    fn aggregate_parent_pct_zero_parent() {
        let e = make_event("x", 0.0, 50, 0);
        let section = ProfileSection::aggregate(&[&e], Some(Duration::ZERO));
        assert_eq!(section.parent_pct, 0.0);
    }

    #[test]
    fn aggregate_parent_pct_none() {
        let e = make_event("x", 0.0, 50, 0);
        let section = ProfileSection::aggregate(&[&e], None);
        assert_eq!(section.parent_pct, 0.0);
    }

    #[test]
    fn profile_section_children_initially_empty() {
        let e = make_event("root", 0.0, 100, 0);
        let section = ProfileSection::aggregate(&[&e], None);
        assert!(section.children.is_empty());
    }
}
