"""
Tests for DNAfold2 core module.
"""

import pytest
from pathlib import Path
import tempfile

from dnafold2.core import DNAFolder, FoldingResult, fold
from dnafold2.config import FoldingConfig

pytestmark = pytest.mark.unit


class TestFoldingResult:
    """Tests for FoldingResult dataclass."""

    def test_result_creation(self):
        """Test creating a FoldingResult."""
        result = FoldingResult(sequence="ATCGATCG")
        assert result.sequence == "ATCGATCG"
        assert result.cg_structures == []
        assert result.all_atom_structures == []
        assert result.elapsed_time == 0.0

    def test_result_repr(self):
        """Test string representation."""
        result = FoldingResult(sequence="ATCGATCGATCGATCGATCGATCG", output_dir=Path("/tmp/test"))
        repr_str = repr(result)
        assert "ATCGATCGATCGATCGATCG..." in repr_str
        assert "n_structures=0" in repr_str


class TestDNAFolder:
    """Tests for DNAFolder class."""

    def test_folder_creation_default_config(self):
        """Test creating folder with default config."""
        # This may fail if binaries aren't installed, but tests the init
        try:
            folder = DNAFolder()
            assert folder.config is not None
            assert folder.config.sampling_method == "remc"
        except FileNotFoundError:
            # Expected if binaries not installed
            pytest.skip("Binaries not installed")

    def test_folder_creation_custom_config(self):
        """Test creating folder with custom config."""
        config = FoldingConfig(sampling_method="sa", folding_steps=100000)
        try:
            folder = DNAFolder(config=config)
            assert folder.config.sampling_method == "sa"
            assert folder.config.folding_steps == 100000
        except FileNotFoundError:
            pytest.skip("Binaries not installed")

    def test_invalid_sequence_rejected(self):
        """Test that invalid sequences are rejected during fold."""
        try:
            folder = DNAFolder()
            with pytest.raises(ValueError, match="Invalid characters"):
                folder.fold("ATCGUXYZ")
        except FileNotFoundError:
            pytest.skip("Binaries not installed")


class TestFoldFunction:
    """Tests for convenience fold function."""

    def test_fold_function_exists(self):
        """Test that fold function is importable."""
        from dnafold2.core import fold

        assert callable(fold)


class TestIntegration:
    """Integration tests (require full installation)."""

    @pytest.mark.smoke
    @pytest.mark.slow
    @pytest.mark.skip(reason="Requires compiled binaries")
    def test_full_folding_pipeline(self):
        """Test complete folding pipeline."""
        config = FoldingConfig(
            sampling_method="remc",
            folding_steps=1000,  # Very short for testing
            optimizing_steps=100,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            result = fold(sequence="ATCCTAGTTATAGGAT", config=config, output_dir=tmpdir)

            assert result.output_dir is not None
            assert result.elapsed_time > 0
