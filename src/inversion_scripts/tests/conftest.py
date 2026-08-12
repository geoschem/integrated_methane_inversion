"""Pytest configuration and shared fixtures for jacobian.py tests."""
import pytest

@pytest.fixture
def temp_dir(tmp_path):
    """Provide a temporary directory for test outputs."""
    return tmp_path
