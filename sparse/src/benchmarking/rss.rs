//! RSS (Resident Set Size) measurement via `/proc/self/status`.

/// Read the process peak RSS (VmHWM) in kilobytes from `/proc/self/status`.
///
/// Returns `None` on non-Linux platforms or if the file cannot be read/parsed.
pub fn read_peak_rss_kb() -> Option<u64> {
    read_proc_status_field("VmHWM:")
}

/// Read the process current RSS (VmRSS) in kilobytes from `/proc/self/status`.
///
/// Returns `None` on non-Linux platforms or if the file cannot be read/parsed.
pub fn read_current_rss_kb() -> Option<u64> {
    read_proc_status_field("VmRSS:")
}

/// Read a kilobyte value from `/proc/self/status` by field prefix.
fn read_proc_status_field(prefix: &str) -> Option<u64> {
    #[cfg(target_os = "linux")]
    {
        let status = std::fs::read_to_string("/proc/self/status").ok()?;
        for line in status.lines() {
            if let Some(rest) = line.strip_prefix(prefix) {
                let trimmed = rest.trim();
                // Format: "12345 kB"
                let kb_str = trimmed.strip_suffix("kB").unwrap_or(trimmed).trim();
                return kb_str.parse::<u64>().ok();
            }
        }
        None
    }
    #[cfg(not(target_os = "linux"))]
    {
        let _ = prefix;
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[cfg(target_os = "linux")]
    fn read_peak_rss_returns_some_on_linux() {
        let rss = read_peak_rss_kb();
        assert!(rss.is_some(), "should read VmHWM on Linux");
        assert!(rss.unwrap() > 0, "peak RSS should be positive");
    }

    #[test]
    #[cfg(not(target_os = "linux"))]
    fn read_peak_rss_returns_none_on_non_linux() {
        let rss = read_peak_rss_kb();
        assert!(rss.is_none());
    }

    #[test]
    #[cfg(target_os = "linux")]
    fn read_current_rss_returns_some_on_linux() {
        let rss = read_current_rss_kb();
        assert!(rss.is_some(), "should read VmRSS on Linux");
        assert!(rss.unwrap() > 0, "current RSS should be positive");
    }

    #[test]
    #[cfg(not(target_os = "linux"))]
    fn read_current_rss_returns_none_on_non_linux() {
        let rss = read_current_rss_kb();
        assert!(rss.is_none());
    }

    #[test]
    #[cfg(target_os = "linux")]
    fn current_rss_lte_peak_rss() {
        let current = read_current_rss_kb().unwrap();
        let peak = read_peak_rss_kb().unwrap();
        assert!(
            current <= peak,
            "current RSS ({current}) should be <= peak RSS ({peak})"
        );
    }
}
