# Project 4 - Ray Tracing
In the previous project you will wrote code to raycast and shade mathematical primitives based on a scene input file into a pixel buffer. In this project you will add recursive raytracing to provide reflection and refraction.
Your program should be resistant to errors and should not segfault or produce undefined behavior. If an error occurs, it should print a message to stderr with “Error:” prefixed to a descriptive error message before returning a non-zero error code. I have a test suite designed to test the robustness of your program.
Your program (raytrace) should have this usage pattern:
raytrace width height input.json output.ppm
