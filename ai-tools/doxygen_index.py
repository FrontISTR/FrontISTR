#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Doxygen Index - Loading and searching Doxygen documentation indices

This module handles loading member and compound indices from Doxygen XML files
and provides search/filter utilities.
"""

from typing import Dict, List, Optional
from pathlib import Path
import xml.etree.ElementTree as ET
import re

# Constants
MIN_RESULTS = 1  # Ensure at least one result if any candidates exist


def _extract_location(element: ET.Element) -> tuple[str, Optional[int]]:
    """Extract location information from XML element."""
    loc = element.find("location")
    file_path = loc.get("file", "") if loc is not None else ""
    line_num = loc.get("line", "") if loc is not None else ""
    try:
        line_int = int(line_num) if line_num else None
    except ValueError:
        line_int = None
    return file_path, line_int


def _extract_brief(element: ET.Element) -> str:
    """Extract brief description from XML element."""
    brief_elem = element.find("briefdescription/para")
    return brief_elem.text.strip() if brief_elem is not None and brief_elem.text else ""


def load_indices(xml_dir: Path) -> tuple[List[Dict], List[Dict]]:
    """
    Load member and compound indices from Doxygen XML files.
    
    Args:
        xml_dir: Path to Doxygen XML output directory
        
    Returns:
        (member_index, compound_index)
    """
    members = []
    compounds = []
    
    for xml_path in xml_dir.rglob("*.xml"):
        if xml_path.name == "index.xml":
            continue
            
        try:
            tree = ET.parse(xml_path)
            root = tree.getroot()
        except ET.ParseError:
            # Skip malformed XML files (may occur in Doxygen output)
            continue
        
        # Process compounds
        for compounddef in root.findall(".//compounddef"):
            compound_name = None  # Reset for each compound
            compound_id = compounddef.get("id", "")
            compound_kind = compounddef.get("kind", "")
            compound_name_elem = compounddef.find("compoundname")
            
            if compound_name_elem is not None and compound_name_elem.text:
                compound_name = compound_name_elem.text.strip()
                
                # Extract location and brief description
                file_path, line_int = _extract_location(compounddef)
                brief = _extract_brief(compounddef)
                
                compounds.append({
                    "compoundname": compound_name,
                    "kind": compound_kind,
                    "file": file_path,
                    "line": line_int,
                    "briefdescription": brief,
                    "compoundid": compound_id,
                    "xml_file": str(xml_path)
                })
            
            # Process members within this compound
            for memberdef in compounddef.findall(".//memberdef"):
                member_id = memberdef.get("id", "")
                member_kind = memberdef.get("kind", "")
                name_elem = memberdef.find("name")
                
                if name_elem is not None and name_elem.text:
                    member_name = name_elem.text.strip()
                    
                    # Extract location and brief description
                    file_path, line_int = _extract_location(memberdef)
                    brief = _extract_brief(memberdef)
                    
                    members.append({
                        "name": member_name,
                        "kind": member_kind,
                        "file": file_path,
                        "line": line_int,
                        # Use compound_name if available (may not be set if compound processing failed)
                        "compoundname": compound_name or "",
                        "briefdescription": brief,
                        "memberid": member_id,
                        "xml_file": str(xml_path)
                    })
    
    return members, compounds


def score_match(item: Dict, tokens: List[str], path_filters: List[str]) -> int:
    """
    Calculate relevance score for search results.
    
    Args:
        item: Item to score (member or compound)
        tokens: Search query tokens
        path_filters: File path filters
        
    Returns:
        Relevance score (higher is better)
    """
    score = 0
    name = (item.get("name") or item.get("compoundname") or "").lower()
    file_path = (item.get("file") or "").lower().replace("\\", "/")
    
    for token in tokens:
        if token in name:
            score += 3
        if token in file_path:
            score += 2
    
    for pattern in path_filters:
        if pattern.lower().replace("\\", "/") in file_path:
            score += 1
    
    return score


def search_and_rank(
    items: List[Dict],
    query: str,
    kind_filter: Optional[List[str]] = None,
    file_pattern: Optional[str] = None,
    limit: int = 50
) -> List[Dict]:
    """
    Common search, filter, and ranking logic.
    
    Args:
        items: List of items to search (members or compounds)
        query: Search query string
        kind_filter: Optional list of kinds to filter by
        file_pattern: Optional file path pattern to filter by
        limit: Maximum number of results
        
    Returns:
        Ranked and filtered list of items
    """
    tokens = [t.lower() for t in re.split(r"\s+", query.strip()) if t.strip()]
    kind_set = set(k.lower() for k in (kind_filter or []))
    
    candidates = []
    for item in items:
        # Filter by kind
        if kind_set and item.get("kind", "").lower() not in kind_set:
            continue
        
        # Filter by file pattern
        if file_pattern:
            file_path = item.get("file", "").lower().replace("\\", "/")
            if file_pattern.lower() not in file_path:
                continue
        
        # Calculate score
        score = score_match(item, tokens, [file_pattern] if file_pattern else [])
        if score > 0 or not tokens:
            candidates.append((score, item))
    
    # Sort by score and name
    name_key = "name" if "name" in (items[0] if items else {}) else "compoundname"
    candidates.sort(key=lambda x: (-x[0], x[1].get(name_key, "")))
    
    # Return top results (remove xml_file from results)
    results = []
    for _, item in candidates[:max(MIN_RESULTS, limit)]:
        result = item.copy()
        result.pop("xml_file", None)
        results.append(result)
    
    return results


def find_members_by_name(member_index: List[Dict], name: str) -> List[Dict]:
    """
    Find members by exact name match (case-insensitive).
    
    Args:
        member_index: Member index to search
        name: Member name to search for
        
    Returns:
        List of matching member records from index
    """
    name_lower = name.lower()
    return [m for m in member_index if m.get("name", "").lower() == name_lower]


def find_compounds_by_name(compound_index: List[Dict], name: str) -> List[Dict]:
    """
    Find compounds by exact name match (case-insensitive).
    
    Args:
        compound_index: Compound index to search
        name: Compound name to search for
        
    Returns:
        List of matching compound records from index
    """
    name_lower = name.lower()
    return [c for c in compound_index if c.get("compoundname", "").lower() == name_lower]
