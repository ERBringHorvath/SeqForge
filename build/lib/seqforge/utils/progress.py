# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Elijah Bring Horvath

import sys
import threading
from time import time
import shutil

class ProgressHandler:
    def __init__(self, total=0, prefix='', mode='bar', width=30, out=sys.stderr):
        self.total = total
        self.prefix = prefix
        self.mode = mode  
        self.requested_width = width
        self.out = out
        self.lock = threading.Lock()
        self.count = 0
        self.start_time = time()

    def _format_time(self, seconds):
        if seconds < 60:
            return f"{seconds:.1f}s"
        m, s = divmod(int(seconds), 60)
        return f"{m}m{s:02d}s"

    def _compute_bar_width(self, extra_len):
        """Compute available bar width given current terminal size and extra characters that will appear."""
        try:
            term_width = shutil.get_terminal_size((80, 20)).columns
        except Exception:
            term_width = 80
        #Reserve space for prefix, percentages, counts, rate, ETA, elapsed, and padding
        reserved = len(self.prefix) + extra_len + 10  # some safety margin
        available = term_width - reserved
        if available < 10:
            return min(self.requested_width, 10)
        return min(self.requested_width, available)

    def update(self, increment=1, verbose_msg=None, current_item=None):
        with self.lock:
            self.count += increment
            if self.mode == 'none':
                return

            elapsed = time() - self.start_time
            rate = self.count / elapsed if elapsed > 0 else 0

            if self.mode == 'verbose':
                if verbose_msg:
                    self.out.write(f"{verbose_msg}\n")
                else:
                    item_info = f" {current_item}" if current_item else ""
                    self.out.write(
                        f"{self.prefix} {self.count}/{self.total}{item_info} "
                        f"({rate:.2f} it/s, elapsed={self._format_time(elapsed)})\n"
                    )
                self.out.flush()
                return

            #bar mode
            if self.total and self.total > 0:
                frac = min(self.count / self.total, 1.0)
                #estimate extra text length aside from bar itself for width calculation
                #e.g., " [ ] 100% (XXX/XXX) 0.00it/s ETA: XX Elapsed: XX"
                extra_desc = f" {int(frac*100):3d}% ({self.count}/{self.total}) {rate:.2f}its/s ETA: {self._format_time(0)}"
                bar_width = self._compute_bar_width(len(extra_desc))
                filled_blocks = int(bar_width * frac)
                empty_blocks = bar_width - filled_blocks
                bar = f"{'█' * filled_blocks}{'-' * empty_blocks}"
                percent = int(frac * 100)
                eta = ((elapsed / frac) - elapsed) if frac > 0 else 0
                eta_str = self._format_time(eta) if frac > 0 else "??"
                line = (
                    f"\r{self.prefix} [{bar}] {percent:3d}% "
                    f"({self.count}/{self.total}) {rate:.2f}it/s "
                )
            else:
                #unknown total: simplified output
                line = (
                    f"\r{self.prefix} {self.count} processed "
                    f"{rate:.2f}it/s"
                )

            self.out.write(line)
            self.out.flush()

            if self.total and self.count >= self.total:
                self.out.flush()

    def finish(self):
        with self.lock:
            if self.mode == 'bar' and self.total and self.count < self.total:
                #force full display
                self.count = self.total
                self.update(0)
