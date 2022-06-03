import os

__version__ = 'Major.Minor.Patch'
version_file = os.path.join(os.path.dirname(__file__), 'VERSION')
if os.path.exists(version_file):
    with open(version_file, 'r') as f:
        __version__ = f.readline().strip()
