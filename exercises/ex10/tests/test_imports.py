import pytest

# Grouping related tests into a class is just for organization purposes
class TestImports:
    def test_magnumopus_import(self):
        """Can the magnumopus package be imported"""
        try:
            import magnumopus
        except:
            pytest.fail(f"Your magnumopus package could not be imported")
    
    def test_read_module_import(self):
        """Can the read module be imported"""
        try:
            from magnumopus import sam
        except:
            pytest.fail(f"Your sam module could not be imported")
    
    def test_read_class_import(self):
        """Can the Read class be imported"""
        try:
            from magnumopus.sam import Read
        except:
            pytest.fail(f"Your Read class could not be imported")
