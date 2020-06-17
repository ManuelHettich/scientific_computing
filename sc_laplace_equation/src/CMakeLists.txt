file(GLOB SOURCES *.cpp)
file(GLOB HEADERS *.h)

add_executable(gradients ${HEADERS} ${SOURCES})

target_link_libraries(gradients glew pmp ${OPENGL_LIBRARIES})

