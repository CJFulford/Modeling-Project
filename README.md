# Boolean Operations On Implicit Models
Created By: Brendan Petras, Conlan Hanwell, and Cody Fulford

This is a tool that implements boolean operations on implicit models.
The purpose behind the project is to assist modeling instructors with teaching their students with an interactive simulation.

## Controls:
    Scroll-Up   Zoom Camera in
    Scroll-Down Zoom Camera out
    LeftClick&Drag  to rotate the camera

    Right Click Select an object
    Right click a second time to select a second object
    Right click on no objects to deselect all objects

    when a single object is selected:
        press g to enter object-move mode
            hold x and drag the mouse across the screens x axis to move the object along the worlds x axis
            hold y and drag the mouse across the screens y axis to move the object along the worlds y axis
            hold z and drag the mouse across the screens y axis to move the object along the worlds y axis
            press g to exit object-move mode
        press r to enter object-rotate mode
            hold x and drag the mouse across the screens x axis to rotate the object about the worlds x axis
            hold y and drag the mouse across the screens y axis to rotate the object about the worlds y axis
            hold z and drag the mouse across the screens y axis to rotate the object about the worlds y axis
            press r to exit object-move mode
        press s to enter object-scale mode
            move the mouse along the screens y axis the scale the object
        press b to break a boolean operation on two objects into its two individual objects
    when two objects are selected
        press U to perform the union operation on the two objects
        press I to perform the intersection operation on the two objects
        press D to perform the subtraction operation on the two objects in the order A-B where A was the first object selected


## The Program   
The program will operate via a real-time ray tracer.

#### Libraries used:
    GLAD V4.5
    GLFW x86
    GLM
    OpenMP

Currently the ray tacer operates only via the CPU, with multithreading via OpenMP.

## Functionality
    Spheres:        Enabled
    Triangles:      Enabled
    Cube:           Enabled
    Cylinder:       Enabled
    Torus:          Enabled

    Operations:     Union
                    Intersction
                    Difference

    Shading:        Blinn-Phong
    Reflections:    Disabled
    Shadows:        Disabled

