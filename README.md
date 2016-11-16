# Project 4 - Ray Tracing

In the previous project you will wrote code to raycast and shade mathematical primitives based on a scene input file into a pixel buffer. 

In this project you will add recursive raytracing to provide reflection and refraction.

Your program should be resistant to errors and should not segfault or produce undefined behavior. If an error occurs, it should print a message to stderr with “Error:” prefixed to a descriptive error message before returning a non-zero error code.

# Example to run the program

From the comand line in the order of :


![alt tag](https://github.com/mhaa54/cs430project04/blob/master/Example%20to%20run%20the%20program.png)

# Specifically, these new properties should be supported for objects:

<br />reflectivity: The amount of reflection (0.0-1.0) <br />
<br />refractivity: The amount of refraction (0.0-1.0) <br />
<br />ior: The index of refraction of the volume<br />

<br />For planes assume that the volume that the index of refraction applies to is on the opposite side of the plane from the camera.<br />

# Scenes for project 4:

![alt tag](https://github.com/mhaa54/cs430project04/blob/master/Scenes.png)
