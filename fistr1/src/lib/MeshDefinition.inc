!=== 1-Dimension mesh ===
  type LINE_DEFINE
!   node
    integer		::	n_n1(2,3)			! Line end node first, then the line midside node
  end type LINE_DEFINE
  type(LINE_DEFINE)		::	LINE
  data LINE%n_n1/	2,3,	1,3,	1,2/

!=== 2-Dimension mesh ===
! Triangle mesh
  type TRI_DEFINE
!   node
	integer   ::	n_n1(4,6)			!	Node_to_adjacent-node
	integer   ::	n_l (2,6)			!	Node_to_line
!   line
	integer		::	l_n (3,3)			!	Line_to_node
  end type TRI_DEFINE
  type(TRI_DEFINE)	::	TRI
  data TRI%n_n1	/	2,3,4,6,		&
                  3,1,5,4,		&
                  1,2,6,5,		&
                  1,2,6,5,		&   ! Mid nodes
                  2,3,4,6,		&
                  3,1,5,4/
  data TRI%n_l	/	1,3,	2,1,	3,2,	&
                  1,0,	2,0,	3,0/
  data TRI%l_n	/	1,2,4,	2,3,5,	3,1,6/

! Quadrilateral mesh
  type QUAD_DEFINE
!   node
	integer		::	n_n1(4,8)			!	Node_to_adjacent-node
	integer		::	n_l (2,8)			!	Node_to_line
!   line
	integer		::	l_n (3,4)			!	Line_to_node
  end type QUAD_DEFINE
  type(QUAD_DEFINE)	::	QUAD
  data QUAD%n_n1/	2,4,	5,8,		&
                  3,1,	6,5,		&
                  4,2,	7,6,		&
                  1,3,	8,7,		&
                  1,2,8,6,      & ! Mid nodes
                  2,3,5,7,      &
                  3,4,6,8,      &
                  4,1,7,5/
  data QUAD%n_l /	1,4, 2,1, 3,2, 4,3,		&
                  1,0, 2,0, 3,0, 4,0/
  data QUAD%l_n /	1,2,5,	2,3,6,	3,4,7,	4,1,8/

!=== 3-Dimension mesh ===
! Tetrahedron mesh
  type TETRA_DEFINE
!   node
	integer		::	n_n1(6,10)			!	Vertex Node
	integer		::	n_l (3,10)			!	Vertex and midside
	integer		::	n_s	(3,10)			!	Node_to_surface
!   line
	integer		::	l_n (3,6)				!	Line_to_node
!	integer,	dimension(4,6)	::	l_4l	!	Line_to_4_line
	integer		::	l_s1(2,6), l_s2(2,6)	!	Line_to_surface
!   surface
	integer		::	s_n (6,4)				!	Surface_to_node
	integer		::	s_l (3,4)				!	Surface_to_line
  end type TETRA_DEFINE
  type(TETRA_DEFINE)	::	TETRA
  data TETRA%n_n1/  2,3,4,	5,7,8,		&
                    3,1,4,	6,5,9,		&
                    1,2,4,	7,6,10,		&
                    1,3,2,	8,10,9,		&
                    1,2,  7,6,8,9,    & ! Mid nodes
                    2,3,  5,7,9,10,   &
                    3,1,	6,5,10,8,   &
                    1,4,	5,9,7,10,   &
                    2,4,	6,10,5,8,   &
                    3,4,	7,8,6,9/
  data TETRA%n_l /  1,3,4,	2,1,5,	3,2,6,	4,6,5,		&
                    1,0,0,	2,0,0,	3,0,0,	4,0,0,	5,0,0,	6,0,0/
  data TETRA%n_s /	1,4,2,	1,2,3,	1,3,4,	4,3,2,		&
                    1,2,0,	1,3,0,	1,4,0,            &
                    2,4,0,	3,2,0,	4,3,0/
!
  data TETRA%l_n /	1,2,5,	2,3,6,	3,1,7,	1,4,8,	2,4,9,	3,4,10/
  data TETRA%l_s1/	1,2,	1,3,	1,4,	2,4,	3,2,	4,3/	! surfaces share the line
  data TETRA%l_s2/	4,3,	2,4,	3,2,	1,3,	1,4,	1,2/	! surfaces connected by the line
!
  data TETRA%s_n /	1,2,3,5,6,7,		&
                    1,4,2,8,9,5,		&
                    2,4,3,9,10,6,		&
                    3,4,1,10,8,7/			! counter-clockwise viewed from inside
  data TETRA%s_l /	1,2,3,	4,5,1,	5,6,2,	6,4,3/

! Hexahedron mesh
  type HEXA_DEFINE
! node
	integer		::	n_n1(6,20), n_n2(3,8)		!	Node_to_node
	integer,	dimension(3,20)	::	n_l				!	Node_to_line
	integer,	dimension(3,26)	::	n_s				!	Node_to_surface
	integer,	dimension(8)    ::	n_dn			!	Diagnal node
	integer,	dimension(8,8)	::	n_o				!   Node order transfer for hexa refinement

! line
	integer,	dimension(3,12)	::	l_n				!	Line_to_node
	integer,	dimension(2,12)	::	l_s1, l_s2		!	Line_to_surface
	integer,	dimension(8,12) ::	l_o				!   Node order transfer for hexa refinement

! surface
	integer,	dimension(9,6)	::	s_n				!	Surface_to_node
	integer,	dimension(4,6)	::	s_l				!	Surface_to_line
	integer,	dimension(8,6)	::	s_o				!   Node order transfer for hexa refinement
  end type HEXA_DEFINE
  type(HEXA_DEFINE)	::	HEXA
! Node
  data HEXA%n_n1/	2,4,5,	9,12,13,	&
                  3,1,6,	10,14,9,	&
                  4,2,7,	11,10,15,	&
                  1,3,8,	12,11,16,	&
                  8,6,1,	20,17,13,	&
                  5,7,2,	17,18,14,	&
                  6,8,3,	18,19,15,	&
                  7,5,4,	19,20,16,	&
                  1,2,	12,10,13,14,	&
                  2,3,	10,11,14,15,	&
                  4,3,	16,15,12,10,	&
                  1,4,	13,16,9,11,		&
                  1,5,	9,17,12,20,		&
                  2,6,	10,18,9,17,		&
                  3,7,	11,19,10,18,	&
                  4,8,	12,20,11,19,	&
                  5,6,	13,14,20,18,	&
                  6,7,	14,15,17,19,	&
                  8,7,	20,18,16,15,	&
                  5,8,	17,19,13,16/
  data HEXA%n_n2/	3,8,6,	4,5,7,	1,6,8,	2,7,5,		&
                  7,2,4,	8,3,1,	5,4,2,	6,1,3/
  data HEXA%n_l /	1,5,9,	8,1,10, 2,8,11,	5,2,12,		&
                  6,4,9,  4,7,10, 7,3,11, 3,6,12,		&
                  1,0,0,	8,0,0,	2,0,0,	5,0,0,		&
                  9,0,0,	10,0,0,	11,0,0,	12,0,0,		&
                  4,0,0,	7,0,0,	3,0,0,	6,0,0/
  data HEXA%n_s/	6,2,4,	6,4,1,	6,1,3,	6,3,2,		&
                  5,4,2,	5,1,4,	5,3,1,	5,2,3,		&
                  6,4,0,	6,1,0,	6,3,0,	6,2,0,		&
                  4,2,0,	1,4,0,	3,1,0,	2,3,0,		&
                  4,5,0,	1,5,0,	5,3,0,	5,2,0,		&
                  6,0,0,	4,0,0,	1,0,0,	3,0,0,	2,0,0,	5,0,0/
  data HEXA%n_dn/	7,8,5,6,3,4,1,2/
! Line
  data HEXA%l_n /	1,2,9,	4,3,11,	8,7,19,	5,6,17,	&
                  1,4,12,	5,8,20,	6,7,18,	2,3,10,	&
                  1,5,13,	2,6,14,	3,7,15,	4,8,11/
  data HEXA%l_s1/	6,4,	3,6,	5,3,	4,5,	&	! surfaces share the line
                  2,6,	5,2,	1,5,	6,1,	&
                  4,2,	1,4,	3,1,	2,3/
  data HEXA%l_s2/	2,1,	2,1,	2,1,	2,1,	&	! surfaces connected by the line
                  4,3,	4,3,	4,3,	4,3,	&
                  6,5,	6,5,	6,5,	6,5/
! Surface
  data HEXA%s_n /	2,6,7,3,	14,18,15,10,	23,	&	! counter-clockwise viewed from inside
                  4,8,5,1,	16,20,13,12,	26,	&
                  3,7,8,4,	15,19,16,11,	24,	&
                  1,5,6,2,	13,17,14,9,		22,	&
                  6,5,8,7,	17,20,19,18,	26,	&
                  1,2,3,4,	9,10,11,12,		21/
  data HEXA%s_l /	10,7,11,8,	12,6,9,5,	&
                  11,3,12,2,	9,4,10,1,	&
                  6,3,7,4,	1,8,2,5/

!  Node order transfer for hexa refinement
    data HEXA%n_o /	6,5,8,7,2,1,4,3,		&
                    7,6,5,8,3,2,1,4,		&
                    8,7,6,5,4,3,2,1,		&
                    5,8,7,6,1,4,3,2,		&
                    4,1,2,3,8,5,6,7,		&
                    1,2,3,4,5,6,7,8,		&
                    2,3,4,1,6,7,8,5,		&
                    3,4,1,2,7,8,5,6/

    data HEXA%l_o /	4,3,7,8,1,2,6,5,		&
                    8,7,6,5,4,3,2,1,		&
                    5,6,2,1,8,7,3,4,		&
                    1,2,3,4,5,6,7,8,		&
                    3,2,6,7,4,1,5,8,		&
                    4,1,2,3,8,5,6,7,		&
                    8,5,1,4,7,6,2,3,		&
                    7,6,5,8,3,2,1,4,		&
                    2,6,7,3,1,5,8,4,		&
                    3,7,8,4,2,6,5,1,		&
                    4,8,5,1,3,7,6,2,		&
                    1,5,6,2,4,8,7,3/

    data HEXA%s_o /	5,1,4,8,6,2,3,7,		&
                    2,6,7,3,1,5,8,4,		&
                    1,5,6,2,4,8,7,3,		&
                    4,3,7,8,1,2,6,5,		&
                    1,2,3,4,5,6,7,8,		&
                    6,5,8,7,2,1,4,3/
! PRISM
  type PRISM_DEFINE
! node
	integer		::	n_n1(6,15)		!	Node_to_node
  integer		::	n_n2(2,6)     !	Vertex node to non-adjacent vertex nodes
	integer,	dimension(3,15)	::	n_l				!	Node_to_line
	integer,	dimension(3,15)	::	n_s				!	Node_to_surface

! line
	integer,	dimension(3,9)	::	l_n				!	Line_to_node
	integer,	dimension(2,9)	::	l_s1, l_s2		!	Line_to_surface

! surface
	integer,	dimension(8,5)	::	s_n				!	Surface_to_node
	integer,	dimension(4,5)	::	s_l				!	Surface_to_line
  end type PRISM_DEFINE
  type(PRISM_DEFINE)	::	PRISM
! node
  data PRISM%n_n1 /	2,3,4,	7,9,10,		&
                    3,1,5,	8,7,11,		&
                    1,2,6,	9,8,12,		&
                    6,5,1,	15,13,10,	&
                    4,6,2,	13,14,11,	&
                    5,4,3,	14,15,12,	&
                    1,2,	9,8,10,11,	&
                    2,3,	7,9,11,12,	&
                    3,1,	8,7,12,10,	&
                    1,4,	7,13,9,15,	&
                    2,5,	8,14,7,13,	&
                    3,6,	9,15,8,14,	&
                    4,5,	10,11,15,14,	&
                    5,6,	11,12,13,15,	&
                    6,4,	12,10,14,13/
  data PRISM%n_n2 /	5,6,	6,4,	4,5,	2,3,	3,1,	1,2/
  data PRISM%n_l  /	1,3,4,	2,1,5,	3,2,6,	9,7,4,	7,8,5,	8,9,6,	&
                    1,0,0,	2,0,0,	3,0,0,	4,0,0,	5,0,0,	6,0,0,	&
                    7,0,0,	8,0,0,	9,0,0/
  data PRISM%n_s  /	5,2,3,	5,3,1,	5,1,2,	4,3,2,	4,1,3,	4,2,1,	&
                    5,3,0,	5,1,0,	5,2,0,	3,2,0,	1,2,0,	2,1,0,	&
                    4,3,0,	4,1,0,	4,2,0/
! line
  data PRISM%l_n  /	1,2,7,	2,3,8,	3,1,9,	1,4,10,	2,5,11,	3,6,12,	&
                    4,5,13,	5,6,14,	6,4,15/
  data PRISM%l_s1 /	5,3,	5,1,	5,2,	3,2,	1,2,	2,1,	&
                    4,3,	4,1,	4,2/
  data PRISM%l_s2 /	2,1,	3,2,	1,2,	4,5,	5,4,	5,4,	&
                    2,1,	3,2,	1,2/
! surface
  data PRISM%s_n  /	2,5,6,3,	11,14,12,8,		&
                    3,6,4,1,	12,15,10,9,		&
                    1,4,5,2,	10,13,11,7,		&
                    4,6,5,		15,14,13,	0,0,&
                    1,2,3,		7,8,9,		0,0/
  data PRISM%s_l  /	5,8,6,2,	6,9,4,3,	4,7,5,1,	&
                    9,8,7,0,	1,2,3,0/
!
! PYRAMID
  type PYRAMID_DEFINE
!   node
	integer		::	n_n1(3,4)		!	Node_to_node 4 nodes on base
  integer   ::  n5_n1(4)    ! On node on top
!
	integer,	dimension(3,4)	::	n_l				!	Node_to_line
	integer,	dimension(4)	::	n5_l				!	Node_to_line
!
	integer,	dimension(3,4)	::	n_s				!	Node_to_surface
	integer,	dimension(4)	::	n5_s			!	Node_to_surface
!   line
	integer,	dimension(2,8)	::	l_n				!	Line_to_node
	integer,	dimension(2,8)	::	l_s     		!	Line_to_surface

!   surface
	integer,	dimension(4,5)	::	s_n				!	Surface_to_node
	integer,	dimension(4,5)	::	s_l				!	Surface_to_line
  end type PYRAMID_DEFINE
  type(PYRAMID_DEFINE)	::	PYRAMID
! node
  data PYRAMID%n_n1 /2,4,5,	3,1,5, 4,2,5, 1,3,5/
  data PYRAMID%n5_n1 /4,3,2,1/

  data PYRAMID%n_l  /1,4,5, 2,1,6, 3,2,7, 4,3,8/
  data PYRAMID%n5_l  /8,7,6,5/

  data PYRAMID%n_s  /5,4,1, 5,1,2, 5,2,3, 5,3,4/
  data PYRAMID%n5_s  /3,2,1,4/

! line
  data PYRAMID%l_n /1,2, 2,3, 3,4, 4,1, 1,5, 2,5, 3,5, 4,5/

  data PYRAMID%l_s /5,1, 5,2, 5,3, 5,4, 1,4, 2,1, 3,2, 4,3/
! surface
  data PYRAMID%s_n /2,1,5,0, 3,2,5,0, 4,3,5,0, 1,4,5,0, 1,2,3,4/

  data PYRAMID%s_l /1,5,6,0, 2,6,7,0, 3,7,8,0, 4,8,5,0, 1,2,3,4/
