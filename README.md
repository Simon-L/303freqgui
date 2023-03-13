# Imgui glfw3 + opengl3 with emscripten support

> Based on imgui example `example_glfw_opengl3`

Moved out of tree and includes workardound for emscripten, currently being investigated

```
commit 773e7e96f95b69605d85e133367ccb4257d6e638 (HEAD -> master)
Author: Simon-L <lumis@xulepth.fr>
Date:   Mon Mar 13 15:24:17 2023 +0100

    Workaround undefined glfw functions in emscripten

diff --git a/backends/imgui_impl_glfw.cpp b/backends/imgui_impl_glfw.cpp
index f0e70a6d..cf6d5386 100644
--- a/backends/imgui_impl_glfw.cpp
+++ b/backends/imgui_impl_glfw.cpp
@@ -536,7 +536,7 @@ static bool ImGui_ImplGlfw_Init(GLFWwindow* window, bool install_callbacks, Glfw
     bd->MouseCursors[ImGuiMouseCursor_NotAllowed] = glfwCreateStandardCursor(GLFW_ARROW_CURSOR);
 #endif
     glfwSetErrorCallback(prev_error_callback);
-#if (GLFW_VERSION_COMBINED >= 3300) // Eat errors (see #5785)
+#if (GLFW_VERSION_COMBINED >= 3300) && !defined(__EMSCRIPTEN__) // Eat errors (see #5785)
     (void)glfwGetError(NULL);
 #endif
 
@@ -669,7 +669,7 @@ static void ImGui_ImplGlfw_UpdateGamepads()
         return;
 
     io.BackendFlags &= ~ImGuiBackendFlags_HasGamepad;
-#if GLFW_HAS_GAMEPAD_API
+#if GLFW_HAS_GAMEPAD_API && !defined(__EMSCRIPTEN__)
     GLFWgamepadstate gamepad;
     if (!glfwGetGamepadState(GLFW_JOYSTICK_1, &gamepad))
         return;
```