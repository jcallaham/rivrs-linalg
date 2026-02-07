//! Peak RSS (Resident Set Size) measurement via `/proc/self/status`.

/// Read the process peak RSS (VmHWM) in kilobytes from `/proc/self/status`.
///
/// Returns `None` on non-Linux platforms or if the file cannot be read/parsed.
pub fn read_peak_rss_kb() -> Option<u64> {
    #[cfg(target_os = "linux")]
    {
        let status = std::fs::read_to_string("/proc/self/status").ok()?;
        for line in status.lines() {
            if let Some(rest) = line.strip_prefix("VmHWM:") {
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
}
