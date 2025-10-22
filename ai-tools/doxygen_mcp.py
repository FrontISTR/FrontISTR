#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Doxygen MCP Server - Fortran/C code documentation access via Model Context Protocol

This server provides tools to query Doxygen-generated XML documentation.
Primary use case: FrontISTR (Fortran/C mixed codebase)

Dependencies: pip install fastmcp
Usage: python doxygen_mcp.py --xml <xml_dir> --html <html_dir>
"""

from __future__ import annotations
import argparse
import sys
import re
from pathlib import Path
import xml.etree.ElementTree as ET
from typing import Dict, List, Optional

from mcp.server.fastmcp import FastMCP, Context

app = FastMCP("doxygen")

# Global data structures
MEMBER_INDEX: List[Dict] = []  # Index of all members
COMPOUND_INDEX: List[Dict] = []  # Index of all compounds
XML_CACHE: Dict[str, ET.Element] = {}  # Cache parsed XML files
XML_DIR: Optional[Path] = None
HTML_ROOT: Optional[Path] = None


def load_indices(xml_dir: Path, html_root: Path) -> tuple[List[Dict], List[Dict]]:
    """
    Load member and compound indices from Doxygen XML files.
    
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
            continue
        
        # Process compounds
        for compounddef in root.findall(".//compounddef"):
            compound_id = compounddef.get("id", "")
            compound_kind = compounddef.get("kind", "")
            compound_name_elem = compounddef.find("compoundname")
            
            if compound_name_elem is not None and compound_name_elem.text:
                compound_name = compound_name_elem.text.strip()
                
                # Location
                loc = compounddef.find("location")
                file_path = loc.get("file", "") if loc is not None else ""
                line_num = loc.get("line", "") if loc is not None else ""
                try:
                    line_int = int(line_num) if line_num else None
                except ValueError:
                    line_int = None
                
                # Brief description
                brief_elem = compounddef.find("briefdescription/para")
                brief = brief_elem.text.strip() if brief_elem is not None and brief_elem.text else ""
                
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
                    
                    # Location
                    loc = memberdef.find("location")
                    file_path = loc.get("file", "") if loc is not None else ""
                    line_num = loc.get("line", "") if loc is not None else ""
                    try:
                        line_int = int(line_num) if line_num else None
                    except ValueError:
                        line_int = None
                    
                    # Brief description
                    brief_elem = memberdef.find("briefdescription/para")
                    brief = brief_elem.text.strip() if brief_elem is not None and brief_elem.text else ""
                    
                    members.append({
                        "name": member_name,
                        "kind": member_kind,
                        "file": file_path,
                        "line": line_int,
                        "compoundname": compound_name if 'compound_name' in locals() else "",
                        "briefdescription": brief,
                        "memberid": member_id,
                        "xml_file": str(xml_path)
                    })
    
    return members, compounds


def score_match(item: Dict, tokens: List[str], path_filters: List[str]) -> int:
    """Calculate relevance score for search results."""
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


@app.tool()
def search_members(
    ctx: Context,
    query: str,
    kind: Optional[List[str]] = None,
    file_pattern: Optional[str] = None,
    limit: int = 50
) -> List[Dict]:
    """
    Search for members (functions, subroutines, variables, types) in Doxygen documentation.
    
    Use this when you need to:
    - Find a function/subroutine by name or keyword
    - Locate variable or type definitions
    - Search for specific procedures in Fortran code
    - Find C functions or global variables
    
    Args:
        query: Search keywords (space-separated). Examples:
               "newton", "matrix stiffness", "damping"
        kind: Filter by member kind. Common values:
              For Fortran: "function", "subroutine", "variable", "typedef"
              For C/C++: "function", "variable", "typedef", "define", "enum"
              Full list: function, variable, typedef, enum, enumvalue, define,
                        signal, slot, property, event, friend, subroutine
        file_pattern: Filter by file path (substring match).
                     Examples: "fistr1/src", "*/dynamic/*", ".f90"
        limit: Maximum number of results (default: 50)
    
    Returns:
        List of basic member information (subset of Doxygen memberdef):
        [
            {
                "name": str,              # Member name (e.g., "fstr_newton")
                "kind": str,              # Member kind (e.g., "function", "subroutine")
                "file": str,              # Source file path
                "line": int,              # Line number in source file
                "compoundname": str,      # Parent compound name (module/namespace, optional)
                "briefdescription": str,  # Brief description (optional)
                "memberid": str           # Doxygen member ID for get_memberdef
            },
            ...
        ]
    
    Examples:
        - Find all subroutines with "newton": kind=["subroutine"], query="newton"
        - Find matrix functions in analysis: query="matrix", file_pattern="*/analysis/*"
    """
    tokens = [t.lower() for t in re.split(r"\s+", query.strip()) if t.strip()]
    kind_set = set(k.lower() for k in (kind or []))
    
    candidates = []
    for member in MEMBER_INDEX:
        # Filter by kind
        if kind_set and member.get("kind", "").lower() not in kind_set:
            continue
        
        # Filter by file pattern
        if file_pattern:
            file_path = member.get("file", "").lower().replace("\\", "/")
            if file_pattern.lower() not in file_path:
                continue
        
        # Calculate score
        score = score_match(member, tokens, [file_pattern] if file_pattern else [])
        if score > 0 or not tokens:
            candidates.append((score, member))
    
    # Sort by score and name
    candidates.sort(key=lambda x: (-x[0], x[1].get("name", "")))
    
    # Return top results
    results = [m for _, m in candidates[:max(1, limit)]]
    
    # Remove xml_file from results (internal use only)
    for r in results:
        r.pop("xml_file", None)
    
    return results


@app.tool()
def search_compounds(
    ctx: Context,
    query: str,
    kind: Optional[List[str]] = None,
    limit: int = 50
) -> List[Dict]:
    """
    Search for compounds (modules, namespaces, types, structs, files) in Doxygen.
    
    Use this when you need to:
    - Find a Fortran module by name
    - Locate type/struct definitions
    - Find C++ namespaces or classes
    - Search for source files
    
    Args:
        query: Search keywords. Examples:
               "step", "tParam", "m_fstr", "solver"
        kind: Filter by compound kind. Common values:
              For Fortran: "namespace" (modules), "struct" (types), "file"
              For C/C++: "namespace", "class", "struct", "file"
              Full list: class, struct, union, interface, namespace,
                        file, dir, page, example, group
        limit: Maximum results (default: 50)
    
    Returns:
        List of basic compound information (subset of Doxygen compounddef):
        [
            {
                "compoundname": str,      # Compound name (e.g., "m_step", "tParamAutoInc")
                "kind": str,              # Compound kind (e.g., "namespace", "struct")
                "file": str,              # Source file (optional)
                "line": int,              # Line number (optional)
                "briefdescription": str,  # Brief description (optional)
                "compoundid": str         # Doxygen compound ID for get_compounddef
            },
            ...
        ]
    
    Examples:
        - Find all Fortran modules: kind=["namespace"]
        - Find type definitions: kind=["struct"], query="tParam"
        - Find files in analysis: kind=["file"], query="analysis"
    """
    tokens = [t.lower() for t in re.split(r"\s+", query.strip()) if t.strip()]
    kind_set = set(k.lower() for k in (kind or []))
    
    candidates = []
    for compound in COMPOUND_INDEX:
        # Filter by kind
        if kind_set and compound.get("kind", "").lower() not in kind_set:
            continue
        
        # Calculate score
        score = score_match(compound, tokens, [])
        if score > 0 or not tokens:
            candidates.append((score, compound))
    
    # Sort by score and name
    candidates.sort(key=lambda x: (-x[0], x[1].get("compoundname", "")))
    
    # Return top results
    results = [c for _, c in candidates[:max(1, limit)]]
    
    # Remove xml_file from results
    for r in results:
        r.pop("xml_file", None)
    
    return results


def extract_text(element: Optional[ET.Element]) -> str:
    """Extract all text content from an XML element."""
    if element is None:
        return ""
    return ET.tostring(element, encoding='unicode', method='text').strip()


def parse_memberdef(memberdef: ET.Element) -> Dict:
    """Parse a memberdef XML element into a dictionary."""
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
        for attr in ["file", "line", "column", "bodyfile", "bodystart", "bodyend"]:
            val = loc.get(attr)
            if val:
                try:
                    loc_dict[attr] = int(val) if attr in ["line", "column", "bodystart", "bodyend"] else val
                except ValueError:
                    loc_dict[attr] = val
        if loc_dict:
            result["location"] = loc_dict
    
    # References (functions this member calls)
    references = []
    for ref in memberdef.findall("references"):
        ref_dict = {"name": extract_text(ref)}
        for attr in ["refid", "compoundref", "startline", "endline"]:
            val = ref.get(attr)
            if val:
                try:
                    ref_dict[attr] = int(val) if attr in ["startline", "endline"] else val
                except ValueError:
                    ref_dict[attr] = val
        if ref_dict.get("name"):
            references.append(ref_dict)
    if references:
        result["references"] = references
    
    # Referenced by (functions that call this member)
    referencedby = []
    for ref in memberdef.findall("referencedby"):
        ref_dict = {"name": extract_text(ref)}
        for attr in ["refid", "compoundref", "startline", "endline"]:
            val = ref.get(attr)
            if val:
                try:
                    ref_dict[attr] = int(val) if attr in ["startline", "endline"] else val
                except ValueError:
                    ref_dict[attr] = val
        if ref_dict.get("name"):
            referencedby.append(ref_dict)
    if referencedby:
        result["referencedby"] = referencedby
    
    return result


@app.tool()
def get_memberdef(
    ctx: Context,
    name: str
) -> Dict:
    """
    Get detailed member information from Doxygen, including call dependencies.
    
    Use this when you need to:
    - Understand what a function/subroutine does
    - See function parameters and return type
    - Find which functions this member calls (call dependencies)
    - Find which functions call this member (caller dependencies)
    - Get documentation and implementation location
    
    This is the primary tool for analyzing function/subroutine dependencies in Fortran/C code.
    
    Args:
        name: Exact member name (case-insensitive).
              Examples: "fstr_newton", "dynamic_mat_ass_load", "NRbound_s"
    
    Returns:
        Dictionary representation of Doxygen <memberdef> element:
        {
            # Basic information
            "id": str,                    # Doxygen member ID
            "kind": str,                  # function|subroutine|variable|typedef|enum|...
            "name": str,                  # Member name
            "type": str,                  # Return type (for functions) or variable type
            "definition": str,            # Full qualified signature
            "argsstring": str,            # Parameter list as string (e.g., "(hecmesh, hecmat, ...)")
            "qualifiedname": str,         # Fully qualified name (e.g., "m_fstr::fstr_newton")
            
            # Documentation
            "briefdescription": str,      # Brief description (optional)
            "detaileddescription": str,   # Detailed description (optional)
            
            # Parameters (for functions/subroutines)
            "param": [                    # Function/subroutine parameters
                {
                    "type": str,          # Parameter type (e.g., "integer(kind=kint)")
                    "declname": str,      # Parameter declaration name
                    "defname": str,       # Parameter definition name
                }
            ],
            
            # Location in source code
            "location": {
                "file": str,              # Source file path
                "line": int,              # Declaration line number
                "bodystart": int,         # Function body start line (optional)
                "bodyend": int            # Function body end line (optional)
            },
            
            # Call dependencies (KEY for dependency analysis)
            "references": [               # Functions/subroutines THIS member calls
                {
                    "refid": str,         # Referenced member ID
                    "name": str,          # Referenced function/subroutine name
                    "compoundref": str,   # Parent compound of referenced member
                    "startline": int,     # Line where call occurs
                    "endline": int
                }
            ],
            "referencedby": [             # Functions/subroutines that call THIS member
                {
                    "refid": str,         # Caller member ID
                    "name": str,          # Caller function/subroutine name
                    "compoundref": str,   # Parent compound of caller
                    "startline": int,
                    "endline": int
                }
            ],
            
            # Multiple match warning (if applicable)
            "_multiple_matches": {        # Present only when multiple members with same name exist
                "total": int,             # Total number of matches found
                "selected": str,          # Description of the selected match
                "others": [str]           # List of other matches (up to 3)
            }
        }
        
        Note: Only elements present in the Doxygen XML will be included.
        For dependency analysis, focus on "references" (calls) and "referencedby" (callers).
        When multiple matches exist, the first match is returned with a warning in _multiple_matches.
    """
    # Find the member in index
    name_lower = name.lower()
    candidates = [m for m in MEMBER_INDEX if m.get("name", "").lower() == name_lower]
    
    if not candidates:
        return {"error": f"Member '{name}' not found"}
    
    # Use first candidate (most common case: only one match)
    member = candidates[0]
    xml_file = member.get("xml_file")
    member_id = member.get("memberid")
    
    if not xml_file or not member_id:
        return {
            "error": "XML details not available",
            "hint": f"Try: get_compounddef(name='{member.get('compoundname', '')}')",
            "basic_info": member
        }
    
    try:
        # Parse XML file
        if xml_file not in XML_CACHE:
            tree = ET.parse(xml_file)
            XML_CACHE[xml_file] = tree.getroot()
        
        root = XML_CACHE[xml_file]
        
        # Find the memberdef by id
        memberdef = root.find(f".//*[@id='{member_id}']")
        if memberdef is None:
            return {"error": "Member not found in XML"}
        
        # Parse and return memberdef details
        result = parse_memberdef(memberdef)
        
        # Add warning if multiple matches exist
        if len(candidates) > 1:
            result["_multiple_matches"] = {
                "total": len(candidates),
                "selected": f"{member.get('compoundname')} in {member.get('file')}",
                "others": [
                    f"{c.get('compoundname')} in {c.get('file')}"
                    for c in candidates[1:4]  # Show up to 3 others
                ]
            }
        
        return result
        
    except Exception as e:
        return {"error": f"Failed to parse XML: {str(e)}"}


def parse_compounddef(compounddef: ET.Element) -> Dict:
    """Parse a compounddef XML element into a dictionary."""
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
        for attr in ["file", "line", "column"]:
            val = loc.get(attr)
            if val:
                try:
                    loc_dict[attr] = int(val) if attr in ["line", "column"] else val
                except ValueError:
                    loc_dict[attr] = val
        if loc_dict:
            result["location"] = loc_dict
    
    # Inner classes/namespaces
    for tag in ["innerclass", "innernamespace", "innerfile", "innerdir"]:
        inners = []
        for inner in compounddef.findall(tag):
            inner_dict = {"name": extract_text(inner)}
            for attr in ["refid", "prot"]:
                val = inner.get(attr)
                if val:
                    inner_dict[attr] = val
            if inner_dict.get("name"):
                inners.append(inner_dict)
        if inners:
            result[tag] = inners
    
    # Base/derived classes
    for tag in ["basecompoundref", "derivedcompoundref"]:
        refs = []
        for ref in compounddef.findall(tag):
            ref_dict = {"name": extract_text(ref)}
            for attr in ["refid", "prot", "virt"]:
                val = ref.get(attr)
                if val:
                    ref_dict[attr] = val
            if ref_dict.get("name"):
                refs.append(ref_dict)
        if refs:
            result[tag] = refs
    
    return result


@app.tool()
def get_compounddef(
    ctx: Context,
    name: str
) -> Dict:
    """
    Get detailed compound (module/namespace/class/file) information from Doxygen.
    
    Use this when you need to:
    - List all functions/subroutines in a Fortran module
    - See all members of a C struct/type
    - Understand the structure of a namespace or class
    - List all symbols defined in a file
    - Find type definitions and their members
    
    For Fortran: compounds are typically modules (namespace kind) or type definitions (struct kind)
    For C/C++: compounds are namespaces, classes, structs, or files
    
    Args:
        name: Compound name. Examples:
              Fortran modules: "m_fstr_nonlinearmethod", "m_step"
              C++ namespaces: "hecmw_util"
              Struct/types: "tParamAutoInc", "fstr_solid"
              Files: "fstr_solve_NonLinear.f90"
    
    Returns:
        Dictionary representation of Doxygen <compounddef> element:
        {
            # Basic information
            "id": str,                    # Doxygen compound ID
            "kind": str,                  # namespace|class|struct|file|group|...
            "language": str,              # Fortran|C|C++|...
            "compoundname": str,          # Compound name
            
            # Documentation
            "briefdescription": str,      # Brief description (optional)
            "detaileddescription": str,   # Detailed description (optional)
            
            # Members organized by sections
            "sectiondef": [               # Member sections
                {
                    "kind": str,          # Section kind (e.g., "public-func", "public-attrib")
                    "memberdef": [        # Members in this section
                        {
                            "id": str,
                            "kind": str,  # function|subroutine|variable|typedef|...
                            "name": str,
                            # ... full memberdef structure
                        }
                    ]
                }
            ],
            
            # Nested compounds
            "innerclass": [               # Nested classes/structs
                {"refid": str, "name": str}
            ],
            "innernamespace": [           # Nested namespaces/modules
                {"refid": str, "name": str}
            ],
            
            # Location
            "location": {
                "file": str,              # Source file path
                "line": int               # Line number (optional)
            }
        }
        
        Note: For type/struct members, look in sectiondef with kind "public-attrib" or similar.
        For module functions, look in sectiondef with kind "func" or "public-func".
    """
    # Find the compound in index
    name_lower = name.lower()
    candidates = [c for c in COMPOUND_INDEX if c.get("compoundname", "").lower() == name_lower]
    
    if not candidates:
        return {"error": f"Compound '{name}' not found"}
    
    compound = candidates[0]
    xml_file = compound.get("xml_file")
    compound_id = compound.get("compoundid")
    
    if not xml_file or not compound_id:
        return {"error": "No XML details available", "basic_info": compound}
    
    try:
        # Parse XML file
        if xml_file not in XML_CACHE:
            tree = ET.parse(xml_file)
            XML_CACHE[xml_file] = tree.getroot()
        
        root = XML_CACHE[xml_file]
        
        # Find the compounddef by id
        compounddef = root.find(f".//*[@id='{compound_id}']")
        if compounddef is None:
            return {"error": "Compound not found in XML", "basic_info": compound}
        
        # Parse and return compounddef details
        return parse_compounddef(compounddef)
        
    except Exception as e:
        return {"error": f"Failed to parse XML: {str(e)}", "basic_info": compound}


def main():
    parser = argparse.ArgumentParser(description="Doxygen MCP Server")
    parser.add_argument("--xml", type=str, default=str(Path.cwd() / "build/doc/xml"),
                       help="Path to Doxygen XML output directory")
    parser.add_argument("--html", type=str, default=str(Path.cwd() / "build/doc/html"),
                       help="Path to Doxygen HTML output directory")
    args = parser.parse_args()
    
    xml_dir = Path(args.xml).resolve()
    html_dir = Path(args.html).resolve()
    
    if not xml_dir.exists():
        print(f"[doxygen-mcp] ERROR: XML directory not found: {xml_dir}", file=sys.stderr)
        sys.exit(2)
    
    if not html_dir.exists():
        print(f"[doxygen-mcp] WARNING: HTML directory not found: {html_dir}", file=sys.stderr)
    
    # Load indices
    global MEMBER_INDEX, COMPOUND_INDEX, XML_DIR, HTML_ROOT
    XML_DIR = xml_dir
    HTML_ROOT = html_dir
    
    print(f"[doxygen-mcp] Loading indices from {xml_dir}...", file=sys.stderr)
    MEMBER_INDEX, COMPOUND_INDEX = load_indices(xml_dir, html_dir)
    
    print(f"[doxygen-mcp] Loaded {len(MEMBER_INDEX)} members, {len(COMPOUND_INDEX)} compounds", file=sys.stderr)
    
    # Run MCP server
    app.run()


if __name__ == "__main__":
    main()
