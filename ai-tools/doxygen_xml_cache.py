#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Doxygen XML Cache - XML file caching and element loading utilities

This module provides XML caching to improve performance when repeatedly
accessing Doxygen XML files.
"""

from typing import Dict, Optional
import xml.etree.ElementTree as ET

# Global XML cache
XML_CACHE: Dict[str, ET.Element] = {}


def load_xml_element(xml_file: str, element_id: str) -> Optional[ET.Element]:
    """
    Load and cache XML file, then find element by ID.
    
    Args:
        xml_file: Path to XML file
        element_id: Doxygen element ID
        
    Returns:
        XML Element or None if not found
    """
    if not xml_file:
        return None
    
    try:
        # Parse and cache XML file
        if xml_file not in XML_CACHE:
            tree = ET.parse(xml_file)
            XML_CACHE[xml_file] = tree.getroot()
        
        root = XML_CACHE[xml_file]
        
        # Find element by id
        return root.find(f".//*[@id='{element_id}']")
        
    except Exception:
        return None


def clear_cache():
    """Clear the XML cache (useful for testing or memory management)."""
    XML_CACHE.clear()
