file(GLOB SOURCES *.cpp)
file(GLOB HEADERS *.h)

add_executable(diffusion ${HEADERS} ${SOURCES})

target_link_libraries(diffusion glew pmp ${OPENGL_LIBRARIES})

