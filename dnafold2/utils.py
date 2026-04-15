"""
DNAfold2 Utility Functions

Helper functions for file handling, process management, and validation.
"""

import logging
import os
import subprocess
import shutil
from pathlib import Path
from typing import Optional, List, Tuple

logger = logging.getLogger(__name__)


def get_package_root() -> Path:
    """Get the root directory of the DNAfold2 package."""
    return Path(__file__).parent.parent


def get_bin_dir() -> Path:
    """Get the directory containing compiled binaries."""
    return get_package_root() / "bin"


def get_src_dir() -> Path:
    """Get the directory containing C source files."""
    return get_package_root() / "src"


def get_data_dir() -> Path:
    """Get the directory containing data files."""
    return get_package_root() / "data"


def check_gcc_available() -> bool:
    """Check if GCC is available for compilation."""
    try:
        result = subprocess.run(
            ["gcc", "--version"],
            capture_output=True,
            text=True
        )
        return result.returncode == 0
    except FileNotFoundError:
        return False


def check_gpp_available() -> bool:
    """Check if G++ is available for compilation."""
    try:
        result = subprocess.run(
            ["g++", "--version"],
            capture_output=True,
            text=True
        )
        return result.returncode == 0
    except FileNotFoundError:
        return False


def compile_source(
    source_file: Path,
    output_file: Path,
    use_openmp: bool = False,
    use_math: bool = True,
    optimize: bool = True
) -> Tuple[bool, str]:
    """Compile a C source file.
    
    Args:
        source_file: Path to .c source file
        output_file: Path for compiled binary
        use_openmp: Enable OpenMP support
        use_math: Link math library
        optimize: Enable -O3 optimization
        
    Returns:
        Tuple of (success, message)
    """
    logger.debug("Compiling %s -> %s", source_file, output_file)
    if not source_file.exists():
        logger.error("Source file not found: %s", source_file)
        return False, f"Source file not found: {source_file}"
    
    # Determine compiler based on file extension
    compiler = "g++" if source_file.suffix == ".cpp" else "gcc"
    
    cmd = [compiler]
    
    if optimize:
        cmd.append("-O3")
    
    cmd.extend(["-Wall", str(source_file), "-o", str(output_file)])
    
    if use_openmp:
        cmd.append("-fopenmp")
    
    if use_math:
        cmd.append("-lm")
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            cwd=source_file.parent
        )
        
        if result.returncode == 0:
            os.chmod(output_file, 0o755)
            logger.info("Successfully compiled %s", source_file.name)
            return True, f"Successfully compiled {source_file.name}"
        else:
            logger.error("Compilation of %s failed: %s", source_file.name, result.stderr[:300])
            return False, f"Compilation failed: {result.stderr}"
    
    except Exception as e:
        logger.error("Compilation error for %s: %s", source_file.name, e)
        return False, f"Compilation error: {str(e)}"


def run_executable(
    executable: Path,
    cwd: Optional[Path] = None,
    timeout: Optional[int] = None,
    env: Optional[dict] = None
) -> Tuple[int, str, str]:
    """Run an executable and capture output.
    
    Args:
        executable: Path to executable
        cwd: Working directory
        timeout: Timeout in seconds
        env: Environment variables
        
    Returns:
        Tuple of (return_code, stdout, stderr)
    """
    logger.debug("Running executable: %s (cwd=%s, timeout=%s)", executable, cwd, timeout)
    if not executable.exists():
        logger.error("Executable not found: %s", executable)
        raise FileNotFoundError(f"Executable not found: {executable}")
    
    if not os.access(executable, os.X_OK):
        os.chmod(executable, 0o755)
        logger.debug("Set execute permission on %s", executable)
    
    run_env = os.environ.copy()
    if env:
        run_env.update(env)
    
    try:
        result = subprocess.run(
            [str(executable)],
            cwd=cwd or executable.parent,
            capture_output=True,
            text=True,
            timeout=timeout,
            env=run_env
        )
        logger.debug("Executable %s returned code %d", executable.name, result.returncode)
        return result.returncode, result.stdout, result.stderr
    
    except subprocess.TimeoutExpired:
        logger.error("Executable %s timed out after %ss", executable.name, timeout)
        return -1, "", "Process timed out"
    except Exception as e:
        logger.error("Executable %s failed: %s", executable.name, e)
        return -1, "", str(e)


def read_sequence_file(seq_file: Path) -> str:
    """Read DNA sequence from a file.
    
    Args:
        seq_file: Path to sequence file
        
    Returns:
        Sequence string
    """
    if not seq_file.exists():
        raise FileNotFoundError(f"Sequence file not found: {seq_file}")
    
    content = seq_file.read_text().strip()
    logger.debug("Read sequence from %s (length=%d)", seq_file, len(content))
    return content


def write_sequence_file(sequence: str, seq_file: Path) -> None:
    """Write DNA sequence to a file.
    
    Args:
        sequence: DNA sequence string
        seq_file: Path to output file
    """
    seq_file.parent.mkdir(parents=True, exist_ok=True)
    seq_file.write_text(sequence.upper().strip() + "\n")
    logger.debug("Wrote sequence to %s", seq_file)


def cleanup_directory(dir_path: Path, keep_results: bool = True) -> None:
    """Clean up a working directory.
    
    Args:
        dir_path: Directory to clean
        keep_results: If True, keep result files
    """
    logger.debug("Cleaning up directory: %s (keep_results=%s)", dir_path, keep_results)
    if not dir_path.exists():
        logger.debug("Directory does not exist, nothing to clean")
        return
    
    if keep_results:
        # Remove only temporary files
        patterns = ["*.o", "a.out", "*.tmp"]
        for pattern in patterns:
            for f in dir_path.glob(pattern):
                f.unlink()
    else:
        shutil.rmtree(dir_path)


def format_pdb_path(base_dir: Path, structure_num: int, structure_type: str = "CG") -> Path:
    """Generate standardized PDB output path.
    
    Args:
        base_dir: Base output directory
        structure_num: Structure number (1-indexed)
        structure_type: "CG" or "All_atom"
        
    Returns:
        Path to PDB file
    """
    subdir = "CG_structure" if structure_type == "CG" else "All_atom_structure"
    filename = f"{structure_type.lower()}_top{structure_num}.pdb"
    return base_dir / subdir / filename
