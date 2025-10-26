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
from pathlib import Path
from typing import Dict, List, Optional

from mcp.server.fastmcp import FastMCP, Context

# Add current directory to path for module imports
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Import internal modules
from doxygen_index import (
    load_indices,
    search_and_rank,
    find_members_by_name,
    find_compounds_by_name
)
from doxygen_parser import parse_memberdef, parse_compounddef
from doxygen_xml_cache import load_xml_element

app = FastMCP("doxygen")

# Global data structures
MEMBER_INDEX: List[Dict] = []  # Index of all members
COMPOUND_INDEX: List[Dict] = []  # Index of all compounds


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
    return search_and_rank(MEMBER_INDEX, query, kind, file_pattern, limit)


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
    return search_and_rank(COMPOUND_INDEX, query, kind, None, limit)


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
    candidates = find_members_by_name(MEMBER_INDEX, name)
    
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
    
    # Load XML element
    memberdef = load_xml_element(xml_file, member_id)
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
    candidates = find_compounds_by_name(COMPOUND_INDEX, name)
    
    if not candidates:
        return {"error": f"Compound '{name}' not found"}
    
    compound = candidates[0]
    xml_file = compound.get("xml_file")
    compound_id = compound.get("compoundid")
    
    if not xml_file or not compound_id:
        return {"error": "No XML details available", "basic_info": compound}
    
    # Load XML element
    compounddef = load_xml_element(xml_file, compound_id)
    if compounddef is None:
        return {"error": "Compound not found in XML", "basic_info": compound}
    
    # Parse and return compounddef details
    return parse_compounddef(compounddef)


def main():
    parser = argparse.ArgumentParser(description="Doxygen MCP Server")
    parser.add_argument("--xml", type=str, default=str(Path.cwd() / "build/doc/xml"),
                       help="Path to Doxygen XML output directory")
    parser.add_argument("--html", type=str, default=str(Path.cwd() / "build/doc/html"),
                       help="Path to Doxygen HTML output directory")
    args = parser.parse_args()
    
    xml_dir = Path(args.xml).resolve()
    html_root = Path(args.html).resolve()
    
    if not xml_dir.exists():
        print(f"[doxygen-mcp] ERROR: XML directory not found: {xml_dir}", file=sys.stderr)
        sys.exit(2)
    
    if not html_root.exists():
        print(f"[doxygen-mcp] WARNING: HTML directory not found: {html_root}", file=sys.stderr)
    
    # Load indices
    global MEMBER_INDEX, COMPOUND_INDEX
    
    print(f"[doxygen-mcp] Loading indices from {xml_dir}...", file=sys.stderr)
    MEMBER_INDEX, COMPOUND_INDEX = load_indices(xml_dir)
    
    print(f"[doxygen-mcp] Loaded {len(MEMBER_INDEX)} members, {len(COMPOUND_INDEX)} compounds", file=sys.stderr)
    
    # Run MCP server
    app.run()


if __name__ == "__main__":
    main()
