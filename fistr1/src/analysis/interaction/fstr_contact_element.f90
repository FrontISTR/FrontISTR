!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief Contact mechanics calculations at element level (single contact pair)
!>
!> Module hierarchy for contact analysis:
!>   Level 1 (Element):  m_fstr_contact_element  - Element calculations for one slave-master pair
!>   Level 2 (Assembly): m_fstr_contact_assembly - Processing all pairs in one tContact object
!>   Level 2 (Search):   m_fstr_contact_search   - Contact detection for all pairs in one tContact
!>   Level 3 (System):   mContact                - Managing all contact objects in the system
!>
!> This module (Level 1) provides:
!>   - Stiffness matrix and force vector computation for a single contact pair
!>   - Tm/Tt matrix computation
!>   - Friction state calculation
!>
!> Module structure:
!>   m_fstr_contact_elem_common - Common utilities (Tm/Tt, friction state, etc.)
!>   m_fstr_contact_elem_slag   - Slag method implementations
!>   m_fstr_contact_elem_alag   - Alag method implementations
!>   m_fstr_contact_element     - Unified interface (this module)
module m_fstr_contact_element
  use m_fstr_contact_elem_common
  use m_fstr_contact_elem_slag
  use m_fstr_contact_elem_alag
  implicit none

  ! Re-export all public functions from sub-modules
  public :: computeTm_Tt
  public :: getTrialFricForceANDcheckFricState
  public :: update_TangentForce
  public :: getContactStiffness_Slag
  public :: getContactNodalForce_Slag
  public :: getTiedStiffness_Slag
  public :: getTiedNodalForce_Slag
  public :: getContactStiffness_Alag
  public :: getContactNodalForce_Alag
  public :: getTiedStiffness_Alag
  public :: getTiedNodalForce_Alag

end module m_fstr_contact_element
