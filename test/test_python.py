import unittest
import os
import shutil
import tempfile
import sys
from pathlib import Path
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
print(f"Parent directory: {parent_dir}")

class TestTemp(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        print(f"\nTEMPDIR: {self.temp_dir}")  # Will show if run correctly
        
        test_path = Path(os.path.join(parent_dir, 'test'))
        self.source_dirs = [item for item in test_path.iterdir() if item.is_dir() and item.name != '__pycache__']
        
        self.copied_dirs = []
        for folder in self.source_dirs:
            dest_folder = os.path.join(self.temp_dir, folder.name)
            shutil.copytree(folder, dest_folder)
            alpha_bridge_path = os.path.join(dest_folder, 'AlphaBridge')
            
            if os.path.exists(alpha_bridge_path):
                shutil.rmtree(alpha_bridge_path)
            
            self.copied_dirs.append(dest_folder)
        
        self.addCleanup(self.cleanup)
            
    def test_example(self):
        self.assertTrue(os.path.exists(self.temp_dir))
        print(self.copied_dirs)

    def cleanup(self):
        if self.temp_dir.exists():
            shutil.rmtree(self.temp_dir)

if __name__ == '__main__':
    unittest.main()