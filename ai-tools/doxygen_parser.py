#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Doxygen XML Parser - Internal utilities for parsing Doxygen XML files

This module contains low-level XML parsing functions used by the Doxygen MCP server.
"""

from typing import Dict, List, Optional, Union
import xml.etree.ElementTree as ET


def extract_text(element: Optional[ET.Element]) -> str:
    """Extract all text content from an XML element."""
    if element is None:
        return ""
    return ET.tostring(element, encoding='unicode', method='text').strip()


def _safe_parse_attr(val: str, attr: str, int_attrs: List[str]) -> Union[int, str]:
    """Parse attribute value to appropriate type (int or str)."""
    if attr in int_attrs:
        try:
            return int(val)
        except ValueError:
            return val
    return val


def _parse_named_refs(parent: ET.Element, tag: str, attrs: List[str], int_attrs: Optional[List[str]] = None) -> List[Dict]:
    """Parse named reference elements with specified attributes."""
    refs = []
    int_attrs = int_attrs or []
    for elem in parent.findall(tag):
        ref_dict = {"name": extract_text(elem)}
        for attr in attrs:
            val = elem.get(attr)
            if val:
                ref_dict[attr] = _safe_parse_attr(val, attr, int_attrs)
        if ref_dict.get("name"):
            refs.append(ref_dict)
    return refs


def parse_memberdef(memberdef: ET.Element) -> Dict:
    """
    Parse a memberdef XML element into a dictionary.
    
    Args:
        memberdef: XML Element representing a Doxygen <memberdef>
        
    Returns:
        Dictionary with parsed member information
    """
    result = {}
    
    # Attributes
    for attr in ["id", "kind", "prot", "static", "const", "explicit", "inline", "virt"]:
        val = memberdef.get(attr)
        if val:
            result[attr] = val
    
    # Simple text elements
    for tag in ["name", "type", "definition", "argsstring", "qualifiedname", "initializer"]:
        elem = memberdef.find(tag)
        if elem is not None:
            text = extract_text(elem)
            if text:
                result[tag] = text
    
    # Description elements
    brief_elem = memberdef.find("briefdescription/para")
    if brief_elem is not None:
        result["briefdescription"] = extract_text(brief_elem)
    
    detailed_elem = memberdef.find("detaileddescription")
    if detailed_elem is not None:
        result["detaileddescription"] = extract_text(detailed_elem)
    
    # Parameters
    params = []
    for param in memberdef.findall("param"):
        param_dict = {}
        for tag in ["type", "declname", "defname", "array", "defval"]:
            elem = param.find(tag)
            if elem is not None:
                text = extract_text(elem)
                if text:
                    param_dict[tag] = text
        if param_dict:
            params.append(param_dict)
    if params:
        result["param"] = params
    
    # Location
    loc = memberdef.find("location")
    if loc is not None:
        loc_dict = {}
        int_attrs = ["line", "column", "bodystart", "bodyend"]
        for attr in ["file", "line", "column", "bodyfile", "bodystart", "bodyend"]:
            val = loc.get(attr)
            if val:
                loc_dict[attr] = _safe_parse_attr(val, attr, int_attrs)
        if loc_dict:
            result["location"] = loc_dict
    
    # References (functions this member calls)
    references = _parse_named_refs(memberdef, "references", 
                                    ["refid", "compoundref", "startline", "endline"],
                                    ["startline", "endline"])
    if references:
        result["references"] = references
    
    # Referenced by (functions that call this member)
    referencedby = _parse_named_refs(memberdef, "referencedby",
                                      ["refid", "compoundref", "startline", "endline"],
                                      ["startline", "endline"])
    if referencedby:
        result["referencedby"] = referencedby
    
    return result


def parse_compounddef(compounddef: ET.Element) -> Dict:
    """
    Parse a compounddef XML element into a dictionary.
    
    Args:
        compounddef: XML Element representing a Doxygen <compounddef>
        
    Returns:
        Dictionary with parsed compound information
    """
    result = {}
    
    # Attributes
    for attr in ["id", "kind", "language", "prot"]:
        val = compounddef.get(attr)
        if val:
            result[attr] = val
    
    # Compound name
    name_elem = compounddef.find("compoundname")
    if name_elem is not None:
        result["compoundname"] = extract_text(name_elem)
    
    # Title (for pages/groups)
    title_elem = compounddef.find("title")
    if title_elem is not None:
        result["title"] = extract_text(title_elem)
    
    # Descriptions
    brief_elem = compounddef.find("briefdescription/para")
    if brief_elem is not None:
        result["briefdescription"] = extract_text(brief_elem)
    
    detailed_elem = compounddef.find("detaileddescription")
    if detailed_elem is not None:
        result["detaileddescription"] = extract_text(detailed_elem)
    
    # Section defs (member groups)
    sections = []
    for sectiondef in compounddef.findall("sectiondef"):
        section_dict = {"kind": sectiondef.get("kind", "")}
        header_elem = sectiondef.find("header")
        if header_elem is not None:
            section_dict["header"] = extract_text(header_elem)
        
        members = []
        for memberdef in sectiondef.findall("memberdef"):
            member_dict = parse_memberdef(memberdef)
            if member_dict:
                members.append(member_dict)
        
        if members:
            section_dict["memberdef"] = members
        sections.append(section_dict)
    
    if sections:
        result["sectiondef"] = sections
    
    # Location
    loc = compounddef.find("location")
    if loc is not None:
        loc_dict = {}
        int_attrs = ["line", "column"]
        for attr in ["file", "line", "column"]:
            val = loc.get(attr)
            if val:
                loc_dict[attr] = _safe_parse_attr(val, attr, int_attrs)
        if loc_dict:
            result["location"] = loc_dict
    
    # Inner classes/namespaces
    for tag in ["innerclass", "innernamespace", "innerfile", "innerdir"]:
        inners = _parse_named_refs(compounddef, tag, ["refid", "prot"])
        if inners:
            result[tag] = inners
    
    # Base/derived classes
    for tag in ["basecompoundref", "derivedcompoundref"]:
        refs = _parse_named_refs(compounddef, tag, ["refid", "prot", "virt"])
        if refs:
            result[tag] = refs
    
    return result
