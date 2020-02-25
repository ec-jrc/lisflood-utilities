import os
import sys


current_dir = os.path.dirname(os.path.abspath(__file__))
src_dir = os.path.join(current_dir, '../src/')
if os.path.exists(src_dir):
    sys.path.append(src_dir)

