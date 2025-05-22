import unittest
import os
import shutil
import tempfile
import sys
from pathlib import Path
import subprocess

parent_dir = Path(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
print(f"Parent directory: {parent_dir}")

class TestTemp(unittest.TestCase):
    def setUp(self):
        
        self.temp_dir = Path(tempfile.mkdtemp())
        test_path = parent_dir / 'test'
        self.source_dirs = [item for item in test_path.iterdir() if item.is_dir() and item.name != '__pycache__']
        
        self.copied_dirs = []
        try:
            for folder in self.source_dirs:
                dest_folder = self.temp_dir / folder.name
                shutil.copytree(folder, dest_folder)
                
                # Remove AlphaBridge if it exists
                alpha_bridge_path = dest_folder / 'AlphaBridge'
                if alpha_bridge_path.exists():
                    shutil.rmtree(alpha_bridge_path)
                
                self.copied_dirs.append(dest_folder)
        except Exception as e:
            self.cleanup()
            raise
        
    def test_example(self):
        """Verify temporary directory and copied folders exist."""
        self.assertTrue(self.temp_dir.exists())
        
        # Additional verification you might want:
        for copied_dir in self.copied_dirs:
            self.assertTrue(Path(copied_dir).exists())
            self.assertTrue(Path(copied_dir).is_dir())
        
    def test_run_script_on_copied_dirs(self):
        
        script_path = parent_dir / "define_interfaces.py"
        script_path = script_path.resolve()
        print(f"Script path: {script_path}") 
        
        self.assertTrue(script_path.exists(), f"Script not found at {script_path}")
        
        
        for input_dir in self.copied_dirs:
            with self.subTest(input_dir=input_dir):  # Creates nested test for each dir
                result = subprocess.run(
                    [
                        sys.executable,
                        str(script_path),
                        "-i",
                        str(input_dir)
                    ],
                    stderr=subprocess.PIPE,
                    stdout=subprocess.PIPE,
                    text=True,
                    env=os.environ,
                    check=True
                )
                
                self._handle_result(result, input_dir)
                
    def _handle_result(self, result, input_dir):
        """Helper to process and validate subprocess results."""
        print(f"\n\033[1mTesting directory: {input_dir}\033[0m")  # Bold text
        
        if result.stderr:
            self._print_error(result.stderr)
        
        # Check for "Finished" in stdout
        self.assertIn("finished", result.stdout, 
                     f"Script did not complete successfully for {input_dir}")
        
        print(f"\033[92m✓ Success: {input_dir}\033[0m")  # Green checkmark
    
    def _print_error(self, error_msg):
        """Print colored error output."""
        red = "\033[91m"
        reset = "\033[0m"
        print(f"{red}═" * 50 + reset)
        print(f"{red}ERROR OUTPUT:{reset}")
        print(f"{red}{error_msg}{reset}")
        print(f"{red}═" * 50 + reset)

    def cleanup(self):
        """Clean up temporary directory."""
        if self.temp_dir.exists():
            try:
                shutil.rmtree(self.temp_dir)
                print(f"Cleaned up: {self.temp_dir}")
            except Exception as e:
                print(f"Warning: Failed to clean up {self.temp_dir}: {e}")

    def tearDown(self):
        """Ensure cleanup runs after each test."""
        self.cleanup()

if __name__ == '__main__':
    unittest.main() 