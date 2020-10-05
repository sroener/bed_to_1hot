"""
Tests the command line interface
"""

from click.testing import CliRunner

from bed_to_1hot import bed_to_1hot
from bed_to_1hot import cli


def test_command_line_interface():
    """Test the CLI."""
    runner = CliRunner()
    #result = runner.invoke(cli.cli, ["-i test_data/encoding_test.bed"])
    #assert result.exit_code == 0
    #assert "bed_to_1hot.cli.cli" in result.output
    help_result = runner.invoke(cli.cli, ["--help"])
    assert help_result.exit_code == 0
    assert "Console script for bed_to_1hot_hd5." in help_result.output
