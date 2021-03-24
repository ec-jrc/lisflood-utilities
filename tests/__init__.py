import os
import sys


current_dir = os.path.dirname(os.path.abspath(__file__))
src_dir = os.path.join(current_dir, '../src/')
if os.path.exists(src_dir):
    sys.path.append(src_dir)

class TestWithCleaner:
    cleanups = []

    def setup_method(self):
        self.cleanups = []

    def teardown_method(self):
        for cleanup_func, args in self.cleanups:
            cleanup_func(*args)
