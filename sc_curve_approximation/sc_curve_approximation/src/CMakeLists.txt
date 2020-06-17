file(GLOB SOURCES *.cpp)
file(GLOB HEADERS *.h)

add_executable(approximation ${HEADERS} ${SOURCES})

target_link_libraries(approximation glew pmp ${OPENGL_LIBRARIES})

