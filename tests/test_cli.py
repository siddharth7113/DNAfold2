"""
Tests for DNAfold2 CLI module.
"""

import pytest
from argparse import Namespace

from dnafold2.cli import create_parser, cmd_validate

pytestmark = pytest.mark.unit


class TestCLI:
    """Tests for command-line interface."""

    def test_parser_creation(self):
        """Test that parser is created successfully."""
        parser = create_parser()
        assert parser is not None

    def test_version_argument(self):
        """Test --version flag."""
        parser = create_parser()
        args = parser.parse_args(["--version"])
        assert args.version is True

    def test_fold_with_sequence(self):
        """Test fold command with inline sequence."""
        parser = create_parser()
        args = parser.parse_args(["fold", "--sequence", "ATCGATCG", "--output", "/tmp/test_output"])
        assert args.command == "fold"
        assert args.sequence == "ATCGATCG"

    def test_fold_with_method(self):
        """Test fold command with method specification."""
        parser = create_parser()
        args = parser.parse_args(["fold", "--sequence", "ATCGATCG", "--method", "sa"])
        assert args.method == "sa"

    def test_fold_with_steps(self):
        """Test fold command with step count."""
        parser = create_parser()
        args = parser.parse_args(["fold", "--sequence", "ATCGATCG", "--steps", "500000"])
        assert args.steps == 500000

    def test_fold_quick_test_flags(self):
        """Test quick test and stage control flags."""
        parser = create_parser()
        args = parser.parse_args(
            [
                "fold",
                "--sequence",
                "ATCGATCG",
                "--quick-test",
                "--stop-after",
                "scoring",
                "--skip-rebuild",
                "--skip-wham",
                "--use-c-stage-tools",
            ]
        )
        assert args.quick_test is True
        assert args.stop_after == "scoring"
        assert args.skip_rebuild is True
        assert args.skip_wham is True
        assert args.use_c_stage_tools is True

    def test_fold_stop_after_initial(self):
        """Test stop-after initial stage choice."""
        parser = create_parser()
        args = parser.parse_args(["fold", "--sequence", "ATCGATCG", "--stop-after", "initial"])
        assert args.stop_after == "initial"

    def test_config_command(self):
        """Test config command."""
        parser = create_parser()
        args = parser.parse_args(
            ["config", "--output", "my_config.dat", "--method", "remc", "--steps", "750000"]
        )
        assert args.command == "config"
        assert args.steps == 750000

    def test_validate_command(self):
        """Test validate command."""
        parser = create_parser()
        args = parser.parse_args(["validate", "ATCGATCG"])
        assert args.command == "validate"
        assert args.sequence == "ATCGATCG"

    def test_info_command(self):
        """Test info command."""
        parser = create_parser()
        args = parser.parse_args(["info"])
        assert args.command == "info"


class TestValidateCommand:
    """Tests for validate command execution."""

    def test_valid_sequence(self):
        """Test validating a valid sequence."""
        result = cmd_validate(Namespace(sequence="ATCGATCG"))
        assert result == 0

    def test_invalid_sequence(self):
        """Test validating an invalid sequence."""
        result = cmd_validate(Namespace(sequence="INVALID"))
        assert result == 1
