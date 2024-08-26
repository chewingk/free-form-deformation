# Free-form Deformation
[Project Report](https://chewingk.github.io/docs/FreeFormDeformation.pdf)

## Introducton
There are many computer graphics techniques developed for computer animation and games. One of the most important techniques has always been object deformation. In both animation and video games, objects would usually have their movements, and it is important to make sure that their shapes do not crack while moving. This could be achieved by artists working on it frame by frame, but that is not cost effective. Researchers have been developing techniques to deform objects, both in 2D and 3D, such as free form deformation (Sederberg & Parry, 1986) which involves fitting the object inside a bounding box, or surface based deformation (Igarashi, Moscovich, & Hughes, 2005).

This project is aiming to implement (Igarashi & Igarashi, 2009), a surface based mesh deformation technique which was established on top of (Igarashi et al., 2005). The program will be written C++, along with the graphics language OpenGL. The input is a 2D triangle mesh of a gingerbread man, which detail is stored in an '.obj' file. The target of this project is to place control points on the mesh, and drag a handle to deform the mesh with the gingerbread man still looking natural.
