# IK_example (0_ikccd)

Practice for some popular IK implemntaion in OpenGL 2.0 

## Installation guide for Windows user (Visual Studio)
1. `git submodule update --init`
2. download glfw [pre-compiled library](https://www.glfw.org/download) and put the uncompressed files under `3rd_party/libglfw`
3. run `cmake -S . -B build` under `/0_ikccd` folder
4. `/0_ikccd/build` will contain the Visual Studio solution files

## Control
1. Using **W S A D** control the horizontal movement of the target 
2. Using **E Q** control the vertical movement of the target
3. Using **ALT + Hold Left Click** control the camera rotation
4. Using **Scroll up/down** control the camera zoom in/out

## TODO
- [x] CCD
    - [x] Basic
    - [x] Hinges
    - [x] Limits
- [ ] FABRIK
- [x] Jacobian Transpose


![title](thumbnail.gif)

# Skeleton IK (1_skeleton_ik)
Applying the CCD IK on rig bones 
## Control
1. Using **W S A D** control the horizontal movement of the target 
2. Using **E Q** control the vertical movement of the target
3. Using **1 2 3 4** switch between four different target IK
4. Using **ALT + Hold Left Click** control the camera
5. Using **Scroll up/down** control the camera zoom in/out
![title](ik_skeleton.gif)

# Jacobian Matrix (using Jacobian Transpose) (2_ik_jacobian)

![jacobian](jacobian.gif)

Practice for some popular IK implemntaion in OpenGL 2.0 
