"""
Tests for DNAfold2 configuration module.
"""

import pytest
from pathlib import Path
import tempfile

from dnafold2.config import FoldingConfig, validate_sequence

pytestmark = pytest.mark.unit


class TestFoldingConfig:
    """Tests for FoldingConfig dataclass."""

    def test_default_config(self):
        """Test default configuration values."""
        config = FoldingConfig()
        assert config.sampling_method == "remc"
        assert config.folding_steps == 150000
        assert config.optimizing_steps == 50000
        assert config.na_concentration == 1000.0
        assert config.mg_concentration == 0.0
        assert config.n_structures == 10

    def test_custom_config(self):
        """Test custom configuration values."""
        config = FoldingConfig(
            sampling_method="sa",
            folding_steps=500000,
            optimizing_steps=100000,
            na_concentration=500.0,
            mg_concentration=10.0,
            n_structures=5,
        )
        assert config.sampling_method == "sa"
        assert config.folding_steps == 500000
        assert config.n_structures == 5

    def test_invalid_sampling_method(self):
        """Test that invalid sampling method raises error."""
        with pytest.raises(ValueError, match="sampling_method must be"):
            FoldingConfig(sampling_method="invalid")

    def test_invalid_negative_steps(self):
        """Test that negative steps raise error."""
        with pytest.raises(ValueError, match="folding_steps must be positive"):
            FoldingConfig(folding_steps=0)

    def test_invalid_negative_concentration(self):
        """Test that negative concentration raises error."""
        with pytest.raises(ValueError, match="na_concentration cannot be negative"):
            FoldingConfig(na_concentration=-100)

    def test_to_dict(self):
        """Test conversion to dictionary."""
        config = FoldingConfig()
        d = config.to_dict()
        assert "sampling_method" in d
        assert "folding_steps" in d
        assert d["sampling_method"] == "remc"

    def test_from_file(self):
        """Test loading configuration from file."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".dat", delete=False) as f:
            f.write("Sampling 1\n")
            f.write("Folding_steps 250000\n")
            f.write("Optimizing_steps 75000\n")
            f.write("CNa 500\n")
            f.write("CMg 5\n")
            f.write("Ncout 8\n")
            f.flush()

            config = FoldingConfig.from_file(f.name)

            assert config.sampling_method == "remc"
            assert config.folding_steps == 250000
            assert config.optimizing_steps == 75000
            assert config.na_concentration == 500.0
            assert config.mg_concentration == 5.0
            assert config.n_structures == 8

            Path(f.name).unlink()

    def test_to_file(self):
        """Test writing configuration to file."""
        config = FoldingConfig(sampling_method="sa", folding_steps=300000)

        with tempfile.NamedTemporaryFile(mode="w", suffix=".dat", delete=False) as f:
            config.to_file(f.name)

            # Read back and verify
            loaded = FoldingConfig.from_file(f.name)
            assert loaded.sampling_method == "sa"
            assert loaded.folding_steps == 300000

            Path(f.name).unlink()


class TestValidateSequence:
    """Tests for sequence validation."""

    def test_valid_sequence(self):
        """Test valid DNA sequence."""
        result = validate_sequence("ATCGATCG")
        assert result == "ATCGATCG"

    def test_lowercase_normalization(self):
        """Test that lowercase is converted to uppercase."""
        result = validate_sequence("atcgatcg")
        assert result == "ATCGATCG"

    def test_mixed_case(self):
        """Test mixed case sequence."""
        result = validate_sequence("AtCgTaCg")
        assert result == "ATCGTACG"

    def test_whitespace_removal(self):
        """Test that whitespace is removed."""
        result = validate_sequence("ATCG ATCG\nTAGC")
        assert result == "ATCGATCGTAGC"

    def test_invalid_characters(self):
        """Test that invalid characters raise error."""
        with pytest.raises(ValueError, match="Invalid characters"):
            validate_sequence("ATCGUATCG")  # U is not valid for DNA

    def test_too_short(self):
        """Test that short sequences raise error."""
        with pytest.raises(ValueError, match="at least 4 nucleotides"):
            validate_sequence("ATG")

    def test_empty_sequence(self):
        """Test that empty sequence raises error."""
        with pytest.raises(ValueError):
            validate_sequence("")
