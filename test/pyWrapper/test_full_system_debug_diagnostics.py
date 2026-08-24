"""Structural contract for full-system debug diagnostics."""

import ast
from pathlib import Path


FULL_SYSTEM_TESTS = Path(__file__).with_name("test_full_system.py")


def test_full_system_tests_define_guarded_local_debug_generators():
    """Require each full-system scenario to retain an opt-in local diagnostic hook."""
    module = ast.parse(FULL_SYSTEM_TESTS.read_text())
    test_functions = [
        node
        for node in module.body
        if isinstance(node, ast.FunctionDef) and node.name.startswith("test_")
    ]

    assert test_functions
    for test_function in test_functions:
        nested_functions = {
            node.name
            for node in test_function.body
            if isinstance(node, ast.FunctionDef)
        }
        calls = [
            node
            for node in ast.walk(test_function)
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
        ]
        call_names = {node.func.id for node in calls}

        assert "generate_debug_data" in nested_functions, test_function.name
        assert "is_debugging" in call_names, test_function.name
        assert "generate_debug_data" in call_names, test_function.name
