# MueLu対応のためのCMakeLists.txt変更例

# 既存のML設定（コメントアウト）
# if(WITH_ML)
#   find_package(ML REQUIRED)
#   target_compile_definitions(hecmw1 PRIVATE HECMW_WITH_ML)
#   target_link_libraries(hecmw1 ${ML_LIBRARIES})
# endif()

# 新しいMueLu設定
option(WITH_MUELU "Enable MueLu (successor to ML)" OFF)

if(WITH_MUELU)
  # Trilinosの検索（MueLu含む）
  find_package(Trilinos REQUIRED COMPONENTS 
    MueLu 
    Tpetra 
    Epetra 
    EpetraExt
    Amesos2 
    Teuchos 
    Ifpack2
    Thyra
  )
  
  # MueLuが利用可能かチェック
  if(NOT ${PROJECT_NAME}_ENABLE_MueLu)
    message(FATAL_ERROR "MueLu was not enabled in this Trilinos build")
  endif()
  
  # Amesos2の利用可能性チェック
  if(${PROJECT_NAME}_ENABLE_Amesos2)
    message(STATUS "Amesos2 found - direct solvers available")
    target_compile_definitions(hecmw1 PRIVATE HAVE_MUELU_AMESOS2)
    
    # 具体的なソルバーの確認
    if(Amesos2_ENABLE_KLU2)
      message(STATUS "KLU2 solver available")
      target_compile_definitions(hecmw1 PRIVATE HAVE_AMESOS2_KLU2)
    endif()
    
    if(Amesos2_ENABLE_MUMPS)
      message(STATUS "MUMPS solver available")
      target_compile_definitions(hecmw1 PRIVATE HAVE_AMESOS2_MUMPS)
    endif()
    
    if(Amesos2_ENABLE_SuperLU)
      message(STATUS "SuperLU solver available") 
      target_compile_definitions(hecmw1 PRIVATE HAVE_AMESOS2_SUPERLU)
    endif()
  endif()
  
  # コンパイル定義
  target_compile_definitions(hecmw1 PRIVATE HECMW_WITH_MUELU)
  
  # インクルードディレクトリ
  target_include_directories(hecmw1 PRIVATE ${Trilinos_INCLUDE_DIRS})
  
  # ライブラリリンク
  target_link_libraries(hecmw1 ${Trilinos_LIBRARIES})
  
  # C++11以上が必要
  target_compile_features(hecmw1 PRIVATE cxx_std_11)
  
  # MueLu wrapper ソースファイルの追加
  target_sources(hecmw1 PRIVATE
    src/solver/precond/ml/hecmw_MueLu_wrapper.cpp
    src/solver/precond/ml/hecmw_MueLu_wrapper_complete.cpp
  )
  
  # MPIが必要
  if(NOT TARGET MPI::MPI_CXX)
    find_package(MPI REQUIRED COMPONENTS CXX)
  endif()
  target_link_libraries(hecmw1 MPI::MPI_CXX)
  
  message(STATUS "MueLu support enabled")
  message(STATUS "Trilinos version: ${Trilinos_VERSION}")
  message(STATUS "Trilinos libraries: ${Trilinos_LIBRARIES}")
  
else()
  message(STATUS "MueLu support disabled")
endif()

# 互換性の確保（MLとMueLuの両方を同時に有効にしない）
if(WITH_ML AND WITH_MUELU)
  message(FATAL_ERROR "Cannot enable both ML and MueLu simultaneously. Please choose one.")
endif()

# デバッグビルド時の追加設定
if(CMAKE_BUILD_TYPE STREQUAL "Debug" AND WITH_MUELU)
  target_compile_definitions(hecmw1 PRIVATE MUELU_DEBUG)
  message(STATUS "MueLu debug mode enabled")
endif()

# インストール設定
if(WITH_MUELU)
  install(FILES 
    src/solver/precond/ml/hecmw_MueLu_wrapper.h
    DESTINATION include
  )
endif()
