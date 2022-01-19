/*
 * Yifei Chen IK practice 2
 * 
 * Control Method 
 * W S A D control the forward/backward/left/right movement of target position 
 * E Q  control the top/down movement of target position
 * 1 2 3 4 switch between different IK target 
 */





#include <iostream>
#include <vector>
#include <algorithm>
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "draw.h"
#include "delfem2/rig_geo3.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/rigv3.h"
#include "delfem2/rig_bvh.h"

namespace dfm2 = delfem2;


dfm2::CQuatd qfv(dfm2::CVec3d& u, dfm2::CVec3d& v)
{
    // a concise version of q from 2 vectors
    u.normalize();
    v.normalize();
    double k = 1.0 + u.dot(v);
    double s = 1 / sqrt(k + k);
    dfm2::CVec3d cross = u.cross(v);
    cross = s * cross;
    return dfm2::CQuatd(k * s, cross.x, cross.y, cross.z);
}


void solveIK(std::vector<dfm2::CRigBone>& aBone, dfm2::CVec3d target, int start_bone, int end_bone)
{
    int num_iter = 3;
    for (unsigned int iter = 0; iter < num_iter; iter++)
    {
        for (unsigned int i = start_bone; i >= end_bone; i--)
        {
            UpdateBoneRotTrans(aBone);
            dfm2::CVec3d jointnow(aBone[i].RootPosition());
            dfm2::CVec3d end(aBone[start_bone + 1].RootPosition());
            dfm2::CVec3d endeff_dir = end - jointnow;
            dfm2::CVec3d target_dir = target - jointnow;

            dfm2::CQuatd r_q;
            r_q = qfv(endeff_dir, target_dir);
            dfm2::CQuatd bone_q(aBone[i].quatRelativeRot[3], aBone[i].quatRelativeRot[0], aBone[i].quatRelativeRot[1], aBone[i].quatRelativeRot[2]);
            dfm2::CQuatd qnew = r_q * bone_q;
            aBone[i].quatRelativeRot[0] = qnew.x;
            aBone[i].quatRelativeRot[1] = qnew.y;
            aBone[i].quatRelativeRot[2] = qnew.z;
            aBone[i].quatRelativeRot[3] = qnew.w;
            UpdateBoneRotTrans(aBone);
        }
    }
}


void draw_coordinate()
{
    ::glDisable(GL_LIGHTING);
    ::glLineWidth(3);
    ::glColor3d(0.1, 0.1, 0.98);
    ::glBegin(GL_LINE_STRIP);
    ::glVertex3d(0, 0, 0);
    ::glVertex3d(0, 0, 20);
    ::glEnd();

    ::glColor3d(0.9, 0.1, 0.11);
    ::glBegin(GL_LINE_STRIP);
    ::glVertex3d(0, 0, 0);
    ::glVertex3d(20, 0, 0);
    ::glEnd();

    ::glColor3d(0.4, 0.4, 0.2);
    ::glBegin(GL_LINE_STRIP);
    ::glVertex3d(0, 0, 0);
    ::glVertex3d(0, 20, 0);
    ::glEnd();
}


void draw_target(dfm2::CVec3d& po, int color)
{
    switch (color)
    {
    case 1 :
        ::glColor3d(0.1, 0.05, 0.08);
        break;
    case 2:
        ::glColor3d(0.9, 0.05, 0.05);
        break;
    case 3:
        ::glColor3d(0.02, 0.05, 0.70);
        break;
    case 4:
        ::glColor3d(0.47, 0.47, 0.02);
        break;
    default:
        ::glColor3d(0.1, 0.05, 0.08);
    }
    ::delfem2::opengl::DrawSphereAt(32, 32, 0.5, po.x, po.y, po.z);
}


static void error_callback(int error, const char* description) {
    fputs(description, stderr);
}

void joystick_callback(int jid, int event)
{
    if (event == GLFW_CONNECTED)
    {
        std::cout << "Controller Connected" << std::endl;
        // The joystick was connected
    }
    else if (event == GLFW_DISCONNECTED)
    {
        std::cout << "Controller Disconnected" << std::endl;
        std::cout << "Recommend connecting a controller for play" << std::endl;

        // The joystick was disconnected
    }
}


int main(int argc, char *argv[]) {
  delfem2::glfw::CViewer3 viewer_source(45);
  viewer_source.width = 1920;
  viewer_source.height = 1080;
  viewer_source.window_title = "Skeleton IK";
  viewer_source.view_rotation = std::make_unique<delfem2::ModelView_Ytop>();

  delfem2::glfw::InitGLOld();

  viewer_source.OpenWindow();
  delfem2::opengl::setSomeLighting();

  glfwSetErrorCallback(error_callback);
  if (!glfwInit())
      exit(EXIT_FAILURE);

  if (!viewer_source.window)
  {
      glfwTerminate();
      exit(EXIT_FAILURE);
  }

  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);

  glfwMakeContextCurrent(viewer_source.window);

  dfm2::CVec3d target{ 10.0,25.0,0.0 };
  dfm2::CVec3d target2{ -5.0,25.0,0.0 };
  dfm2::CVec3d target3{ 7.0,0.0,0.0 };
  dfm2::CVec3d target4{ -2.0,0.0,0.0 };
  std::vector<dfm2::CVec3d> target_group{ target, target2, target3, target4 };

  double ws_forward = 0.0f;
  double ad_leftright = 0.0f;
  double qe_topdown = 0.0f;
  std::string path_bvh = std::string(PATH_SOURCE_DIR) + "/../data/LocomotionFlat01_000.bvh";

  std::vector<dfm2::CRigBone> aBone;
  std::vector<dfm2::CChannel_BioVisionHierarchy> aChannelRotTransBone;
  std::vector<double> aValRotTransBone;

  size_t nframe = 0;
  double frame_time;
  std::string header_bvh;
  std::vector<double> vec_bvh_time_series_data;
  Read_BioVisionHierarchy(
      aBone, 
      aChannelRotTransBone, 
      nframe, 
      frame_time, 
      aValRotTransBone, 
      header_bvh, 
      path_bvh);

  SetPose_BioVisionHierarchy(
      aBone, aChannelRotTransBone,
      aValRotTransBone.data());

  UpdateBoneRotTrans(aBone);
  Floor floor{ 100, -0.1, };
  int target_index = 0;
  int hightlight_bone = 0;
  // there are 4 ik point can be selected 
  
  dfm2::CVec2i left_hand{22,20};
  dfm2::CVec2i right_hand{31,29};
  dfm2::CVec2i left_foot{4,2};
  dfm2::CVec2i right_foot{10,8};
  std::vector<dfm2::CVec2i> ik_group{ left_hand ,right_hand,left_foot,right_foot };


  while (!glfwWindowShouldClose(viewer_source.window))
  {
      glfwMakeContextCurrent(viewer_source.window);
      viewer_source.DrawBegin_oldGL();

      
      int present = 0;
      bool check_controller = false;
      float drop_speed = 0.2;
      float gamepad_x = 0.0f;
      float gamepad_y = 0.0f;
      float gamepad2_x = 0.0f;
      float gamepad2_y = 0.0f;

      glfwSetJoystickCallback(joystick_callback);
      if (glfwJoystickPresent(GLFW_JOYSTICK_1) == GLFW_TRUE)
      {

          present = glfwJoystickPresent(GLFW_JOYSTICK_1);

          if (present == 1)
          {
              check_controller = true;
              int axesCount;
              const float* axes = glfwGetJoystickAxes(GLFW_JOYSTICK_1, &axesCount);
              gamepad_x = axes[0];
              gamepad_y = axes[1];
              gamepad2_x = abs(axes[2]) < 0.05 ? 0.0f : axes[2] / 100;
              gamepad2_y = abs(axes[3]) < 0.05 ? 0.0f : axes[3] / 100;
              float gamepad_mag = sqrtf(gamepad_x * gamepad_x + gamepad_y * gamepad_y);

              int triggerRT = axes[5];
              int triggerLT = axes[4];
              if (gamepad_mag > 0.2f)
              {
                  float gamepaddirx = gamepad_x / gamepad_mag;
                  float gamepaddiry = gamepad_y / gamepad_mag;

                  target_group[target_index].x += gamepaddirx;
                  target_group[target_index].z += gamepaddiry;

              }
              if (triggerLT == 1) {
                  target_group[target_index].y += -drop_speed;
              }
              if (triggerRT == 1) {
                  target_group[target_index].y += drop_speed;
              }
          }
      }
      else
      {
          int state_1 = glfwGetKey(viewer_source.window, GLFW_KEY_1);
          int state_2 = glfwGetKey(viewer_source.window, GLFW_KEY_2);
          int state_3 = glfwGetKey(viewer_source.window, GLFW_KEY_3);
          int state_4 = glfwGetKey(viewer_source.window, GLFW_KEY_4);

          int state_F = glfwGetKey(viewer_source.window, GLFW_KEY_F);
          int state_w = glfwGetKey(viewer_source.window, GLFW_KEY_W);
          int state_s = glfwGetKey(viewer_source.window, GLFW_KEY_S);
          int state_a = glfwGetKey(viewer_source.window, GLFW_KEY_A);
          int state_d = glfwGetKey(viewer_source.window, GLFW_KEY_D);
          int state_q = glfwGetKey(viewer_source.window, GLFW_KEY_Q);
          int state_e = glfwGetKey(viewer_source.window, GLFW_KEY_E);

          if (state_F == GLFW_PRESS)
              hightlight_bone += 1;

          if (state_1 == GLFW_PRESS)
              target_index = 0;
          if (state_2 == GLFW_PRESS)
              target_index = 1;
          if (state_3 == GLFW_PRESS)
              target_index = 2;
          if (state_4 == GLFW_PRESS)
              target_index = 3;

          if (state_w == GLFW_PRESS)
          {
              ws_forward = -0.3f;
          }
          if (state_s == GLFW_PRESS)
          {
              ws_forward = 0.3f;
          }
          if ((state_s == GLFW_RELEASE) && (state_w == GLFW_RELEASE))
          {
              ws_forward = 0.0f;
          }
          if (state_a == GLFW_PRESS)
          {
              ad_leftright = -0.3f;
          }
          if (state_d == GLFW_PRESS)
          {
              ad_leftright = 0.3f;
          }
          if ((state_d == GLFW_RELEASE) && (state_a == GLFW_RELEASE))
          {
              ad_leftright = 0.0f;
          }
          if (state_q == GLFW_PRESS)
          {
              qe_topdown = 0.3f;
          }
          if (state_e == GLFW_PRESS)
          {
              qe_topdown = -0.3f;
          }
          if ((state_e == GLFW_RELEASE) && (state_q == GLFW_RELEASE))
          {
              qe_topdown = 0.0f;
          }
      }
      target_group[target_index].x += ad_leftright;
      target_group[target_index].y += qe_topdown;
      target_group[target_index].z += ws_forward;


      SetPose_BioVisionHierarchy(
          aBone, aChannelRotTransBone,
          aValRotTransBone.data());
      solveIK(aBone, target_group[target_index], ik_group[target_index].x, ik_group[target_index].y);
      for (unsigned int i = 0; i < 4; i++)
      {
          if (i != target_index)
          {
              solveIK(aBone, target_group[i], ik_group[i].x, ik_group[i].y);
          }
      }

      ::glLineWidth(1);
      dfm2::opengl::DrawBone_Octahedron(
          aBone,
          ik_group[target_index].x, -1,
          0.1, 1.0);
      DrawFloorShadow(
          aBone, floor,
          -0.1);
      { // draw floor (stencil 1)
          glEnable(GL_STENCIL_TEST);
          glStencilFunc(GL_ALWAYS, 1, static_cast<GLuint>(~0));
          glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);
          floor.draw_checkerboard();
      }
      draw_target(target_group[0], 1);
      draw_target(target_group[1], 2);
      draw_target(target_group[2], 3);
      draw_target(target_group[3], 4);
      draw_coordinate();
      //std::cout << hightlight_bone << std::endl;
      viewer_source.SwapBuffers();
      //glfwSwapBuffers(viewer_source.window);
      glfwPollEvents();
  }

  glfwDestroyWindow(viewer_source.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
