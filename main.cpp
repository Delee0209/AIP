#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#define STB_IMAGE_RESIZE_IMPLEMENTATION

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>
#include <string.h>
#include <Windows.h>
#include <time.h>
#include <vector>

// GLEW - for openGL's functionality
#include <GL/glew.h>

// GLFW - for drawing Window in screen
#include <GLFW/glfw3.h>

// Imgui - tool for ui
#include "imgui.h"
#include "imgui_internal.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

// stb - for load image
#include "stb_image.h"
#include "stb_image_write.h"
#include "stb_image_resize.h"

// NFDE - for file explorer
#include "nfd.h"

// PPMIO - for write out PPM
#include "ppm_io.h"

# define PI 3.14159265358979323846

GLFWwindow *Window;
int Window_x = 800, Window_y = 500;

static void HelpMarker(const char* desc)
{
    ImGui::TextDisabled("(?)");
    if (ImGui::IsItemHovered())
    {
        ImGui::BeginTooltip();
        ImGui::PushTextWrapPos(ImGui::GetFontSize() * 35.0f);
        ImGui::TextUnformatted(desc);
        ImGui::PopTextWrapPos();
        ImGui::EndTooltip();
    }
}

GLubyte* Load_Texture(int image_width, int image_height, GLuint &image_texture)
{
    GLuint fbo;
    glGenFramebuffers(1, &fbo);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, image_texture, 0);

    int data_size = image_width * image_width * 4;
    GLubyte* pixels = new GLubyte[image_width * image_height * 4];
    glReadPixels(0, 0, image_width, image_height, GL_RGBA, GL_UNSIGNED_BYTE, pixels);

    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glDeleteFramebuffers(1, &fbo);
    return pixels;
}

void Save_Texture(int image_width, int image_height, GLuint &image_texture, GLubyte* pixels)
{
    // Create a OpenGL texture identifier
    glDeleteTextures(1, &image_texture);
    glGenTextures(1, &image_texture);
    glBindTexture(GL_TEXTURE_2D, image_texture);

    // Setup filtering parameters for input image
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE); // This is required on WebGL for non power-of-two textures
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE); // Same

    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, image_width, image_height, 0, GL_RGBA, GL_UNSIGNED_BYTE, pixels);
}

bool LoadImage(std::string &InputFileName, int &input_image_width, int &input_image_height, int &output_image_width, int &output_image_height, GLuint &input_image_texture, GLuint &output_image_texture)
{
    unsigned char* input_image_data = stbi_load(InputFileName.c_str(), &input_image_width, &input_image_height, NULL, 4);
    if (input_image_data == NULL)
    {
        printf("fail to load default image\n");
        return false;
    }

    Save_Texture(input_image_width,input_image_height,input_image_texture,input_image_data);

    output_image_height = input_image_height;
    output_image_width = input_image_width;

    Save_Texture(output_image_width,output_image_height,output_image_texture,input_image_data);

    stbi_image_free(input_image_data);
    return true;
}

void SaveImage(std::string &File_Path, int &output_image_width, int &output_image_height, GLuint &output_image_texture)
{
    GLubyte* pixels = Load_Texture(output_image_width, output_image_height, output_image_texture);

    char format[5] = {'\0'};
    strncpy(format, File_Path.c_str()+File_Path.size()-4,4);
    std::string format_comp = std::string(format);

    if(format_comp == ".jpg")
    {
        //printf("jpg");
        stbi_write_jpg( File_Path.c_str(), output_image_width, output_image_height, 4, pixels, 0);
    }
    else if(format_comp == ".png")
    {
        //printf("png");
        stbi_write_png( File_Path.c_str(), output_image_width, output_image_height, 4, pixels, 0);
    }
    else if(format_comp == ".bmp")
    {
        //printf("bmp");
        stbi_write_bmp( File_Path.c_str(), output_image_width, output_image_height, 4, pixels );
    }
    else if(format_comp == ".ppm")
    {
        //printf("ppm");
        PPM ppm(pixels,output_image_height,output_image_width,255,"P6");
        ppm.write(File_Path.c_str());
    }
    free(pixels);
}

nfdfilteritem_t filterItem[5] = { { "All Support Image File", "jpg,png,ppm,bmp,jpeg" }, { "Image File JPEG", "jpg,jpeg" },{ "Image File PNG", "png" },{ "Image File PPM", "ppm" },{ "Image File BMP", "bmp" }};

bool BrowseFile(std::string &File_Path)
{
    NFD_Init();

    nfdchar_t *outPath;
    nfdresult_t result = NFD_OpenDialog(&outPath, filterItem, 5, NULL);
    if ( result == NFD_OKAY )
    {
        //puts("Success!");
        //puts(outPath);
        File_Path = std::string(outPath);
        NFD_FreePath(outPath);
        return true;
    }
    else if ( result == NFD_CANCEL )
    {
        //puts("User pressed cancel.");
    }
    else
    {
        printf("Error: %s\n", NFD_GetError() );
    }

    NFD_Quit();
    return false;
}

bool BrowseFileSave(std::string &File_Path, int output_image_width, int output_image_height, GLint output_image_texture)
{
    NFD_Init();

    nfdchar_t *outPath;
    nfdresult_t result = NFD_SaveDialog(&outPath, filterItem+1, 4, NULL, "default");
    if ( result == NFD_OKAY )
    {
        //puts("Success!");
        //puts(outPath);
        File_Path = std::string(outPath);

        NFD_FreePath(outPath);
        return true;
    }
    else if ( result == NFD_CANCEL )
    {
        //puts("User pressed cancel.");
    }
    else
    {
        printf("Error: %s\n", NFD_GetError() );
    }

    NFD_Quit();
    return false;
}

void AverageGrayScale(int &image_width, int &image_height, GLubyte* pixels)
{
    for(int i = 0;i < image_height;i ++)
    {
        for(int j = 0;j < image_width;j ++)
        {
            float grayscale = 0;
            for(int k = 0;k < 3;k ++)
            {
                grayscale += pixels[4*(i * image_width + j) + k];
            }
            grayscale /= 3.f;
            for(int k = 0;k < 3;k ++)
            {
                pixels[4*(i * image_width + j) + k] = grayscale * (pixels[4*(i * image_width + j) + 3]/255);
            }
        }
    }
}

void WeightedGrayScale(int &image_width, int &image_height, GLubyte* pixels)
{
    for(int i = 0;i < image_height;i ++)
    {
        for(int j = 0;j < image_width;j ++)
        {
            float grayscale = 0;
            grayscale += 0.299*pixels[4*(i * image_width + j) + 0];
            grayscale += 0.587*pixels[4*(i * image_width + j) + 1];
            grayscale += 0.114*pixels[4*(i * image_width + j) + 2];
            for(int k = 0;k < 3;k ++)
            {
                pixels[4*(i * image_width + j) + k] = grayscale * (pixels[4*(i * image_width + j) + 3]/255);
            }
        }
    }
}

float GaussianDistribution[100];
int GaussianDistributionBinCount = 100;
std::vector<float> z;

void WhiteNoise(int &image_width, int &image_height, GLubyte* pixels, float StandardDeviation)
{
    z.clear();
    for(int i = 0;i < image_height;i ++)
    {
        for(int j = 0;j < image_width;j ++)
        {
            if(j % 2 == 0)
            {
                float r = rand() / (RAND_MAX + 1.0), phi = rand() / (RAND_MAX + 1.0);
                float z1 = StandardDeviation*cos(2*PI*phi)*sqrt(-2*log(r)),z2=StandardDeviation*sin(2*PI*phi)*sqrt(-2*log(r));
                for(int k = 0;k < 3;k++)
                {
                    int tmp = pixels[4*(i * image_width + j) + k];
                    tmp += z1*255;
                    if(tmp > 255) pixels[4*(i * image_width + j) + k] = 255;
                    else if(tmp < 0) pixels[4*(i * image_width + j) + k] = 0;
                    else pixels[4*(i * image_width + j) + k] = tmp;
                    tmp = pixels[4*(i * image_width + j+1) + k];
                    tmp += z2*255;
                    if(tmp > 255) pixels[4*(i * image_width + j+1) + k] = 255;
                    else if(tmp < 0) pixels[4*(i * image_width + j+1) + k] = 0;
                    else pixels[4*(i * image_width + j+1) + k] = tmp;
                }
                z.emplace_back(z1);
                z.emplace_back(z2);
            }
        }
        if(image_width%2==1)
        {
            int j = image_width -1;
            float r = rand() / (RAND_MAX + 1.0), phi = rand() / (RAND_MAX + 1.0);
            float z1 = StandardDeviation*cos(2*PI*phi)*sqrt(-2*log(r));
            for(int k = 0;k < 3;k++)
            {
                int tmp = pixels[4*(i * image_width + j) + k];
                tmp += z1*255;
                if(tmp > 255) pixels[4*(i * image_width + j) + k] = 255;
                else if(tmp < 0) pixels[4*(i * image_width + j) + k] = 0;
                else pixels[4*(i * image_width + j) + k] = tmp;
            }
            z.emplace_back(z1);
        }
    }
}

void WhiteNoiseRGB(int &image_width, int &image_height, GLubyte* pixels, float StandardDeviation)
{
    z.clear();
    for(int i = 0;i < image_height;i ++)
    {
        for(int j = 0;j < image_width;j ++)
        {
            for(int k = 0;k < 3;k++)
            {
                if(j % 2 == 0)
                {
                    float r = rand() / (RAND_MAX + 1.0), phi = rand() / (RAND_MAX + 1.0);
                    float z1 = StandardDeviation*cos(2*PI*phi)*sqrt(-2*log(r)),z2=StandardDeviation*sin(2*PI*phi)*sqrt(-2*log(r));

                    int tmp = pixels[4*(i * image_width + j) + k];
                    tmp += z1*255;
                    if(tmp > 255) pixels[4*(i * image_width + j) + k] = 255;
                    else if(tmp < 0) pixels[4*(i * image_width + j) + k] = 0;
                    else pixels[4*(i * image_width + j) + k] = tmp;
                    tmp = pixels[4*(i * image_width + j+1) + k];
                    tmp += z2*255;
                    if(tmp > 255) pixels[4*(i * image_width + j+1) + k] = 255;
                    else if(tmp < 0) pixels[4*(i * image_width + j+1) + k] = 0;
                    else pixels[4*(i * image_width + j+1) + k] = tmp;

                    z.emplace_back(z1);
                    z.emplace_back(z2);
                }
            }
        }
        if(image_width%2==1)
        {
            int j = image_width -1;
            for(int k = 0;k < 3;k++)
            {
                float r = rand() / (RAND_MAX + 1.0), phi = rand() / (RAND_MAX + 1.0);
                float z1 = StandardDeviation*cos(2*PI*phi)*sqrt(-2*log(r));

                int tmp = pixels[4*(i * image_width + j) + k];
                tmp += z1*255;
                if(tmp > 255) pixels[4*(i * image_width + j) + k] = 255;
                else if(tmp < 0) pixels[4*(i * image_width + j) + k] = 0;
                else pixels[4*(i * image_width + j) + k] = tmp;

                z.emplace_back(z1);
            }
        }
    }
}

void HaarWavelet(int &image_width, int &image_height, GLubyte* pixels, int level, float scalar)
{
    int processRange_width = image_width,processRange_height = image_height;
    bool odd_flag; // 1->odd 0->even
    for(int k = 0;k < level;k ++)
    {
        // 1D Haar Wavelet transform along x
        odd_flag = processRange_width%2==1;
        for(int i = 0;i < processRange_height;i++)
        {
            for(int j = 0;j < processRange_width/2;j ++)
            {
                float S = 0,D = 0;
                // low pass
                S += (float)pixels[4*(i * image_width + 2 * j)];
                S += (float)pixels[4*(i * image_width + 2 * j + 1)];
                S /= 2;
                // high pass
                D += (float)pixels[4*(i * image_width + 2 * j)];
                D -= (float)pixels[4*(i * image_width + 2 * j + 1)];
                D = abs(D);
                D /= 2; // ?
                for(int e = 0;e < 3;e ++)pixels[4*(i * image_width + j)+e] = S;
                pixels[4*(i * image_width + j)+3] = D * scalar;
            }
            for(int j = 0;j < processRange_width/2;j ++)
            {
                for(int e = 0;e < 3;e ++)pixels[4*(i * image_width + j + (processRange_width+1)/2)+e] = pixels[4*(i * image_width + j)+3];
                pixels[4*(i * image_width + j)+3] = 255;
                pixels[4*(i * image_width + j + (processRange_width+1)/2)+3]=255;
            }
            if(odd_flag)
            {
                for(int e = 0;e < 3;e++)
                {
                    pixels[4*(i * image_width + processRange_width/2)+e] = pixels[4*(i * image_width + processRange_width - 1)+e];
                }
                pixels[4*(i * image_width + processRange_width/2)+3]=255;
            }
        }
        // 1D Haar Wavelet transform along y
        odd_flag = processRange_height%2==1;
        for(int i = 0;i < processRange_width;i++)
        {
            for(int j = 0;j < processRange_height/2;j++)
            {
                float S = 0,D = 0;
                // low pass
                S += (float)pixels[4*(2*j * image_width + i)];
                S += (float)pixels[4*((2*j+1) * image_width + i)];
                S /= 2;
                // high pass
                D += (float)pixels[4*(2*j * image_width + i)];
                D -= (float)pixels[4*((2*j+1) * image_width + i)];
                D = abs(D);
                D /= 2; // ?

                for(int e = 0;e < 3;e ++)pixels[4*(j * image_width + i)+e] = S;
                pixels[4*(j * image_width + i)+3] = D * scalar;
            }
            for(int j = 0;j < processRange_height/2;j ++)
            {
                for(int e = 0;e < 3;e ++) pixels[4 * ((j + (processRange_height + 1) / 2) * image_width + i) + e] = pixels[4*(j * image_width + i)+3];
                pixels[4*(j * image_width + i)+3] = 255;
                pixels[4 * ((j + (processRange_height + 1) / 2) * image_width + i)+3]=255;
            }
            if(odd_flag)
            {
                for(int e = 0;e < 3;e++)
                {
                    pixels[4*((processRange_height/2)+e * image_width + i) + e] = pixels[4*((processRange_height - 1) * image_width + i)+e];
                }
            }
        }
        processRange_width = (processRange_width+1)/2;
        processRange_height = (processRange_height+1)/2;
    }
}

void HistogramEqualization(int &image_width, int &image_height, GLubyte* pixels, int mode) // mode 0 -> grey, 1 -> RGB
{
    for(int k = 0;k < 3;k ++)
    {
        // Calculate Histogram
        std::vector<int> Histogram(256, 0);
        for(int i = 0 ;i < image_height;i ++)
        {
            for(int j = 0;j < image_width;j ++)
            {
                Histogram[pixels[4 * (i * image_width + j) + k]] ++;
            }
        }
        // Calculate CDF of Histogram
        std::vector<int> Histogram_CDF(256, 0);
        int G_min = -1;
        Histogram_CDF[0] = Histogram[0];
        for(int i = 1;i < 256;i ++)
        {
            Histogram_CDF[i] = Histogram_CDF[i - 1] + Histogram[i];
            if(G_min == -1 && Histogram[i] > 0) G_min = i;
        }
        int Histogram_min = Histogram_CDF[G_min];
        // Calculate Histogram Equalization
        std::vector<int> Histogram_result(256, 0);
        for(int i = 0;i < 256;i ++)
        {
            Histogram_result[i] = round(((double)((double)Histogram_CDF[i] - (double)Histogram_min)/((double)image_width * (double)image_height - (double)Histogram_min)) * 255.f);
            if(Histogram_result[i] < 0) Histogram_result[i] = 0;
            if(Histogram_result[i] > 255) Histogram_result[i] = 255;
        }
        // write result
        if(mode == 0)
        {
            for(int i = 0 ;i < image_height;i ++)
            {
                for(int j = 0;j < image_width;j ++)
                {
                    int result = Histogram_result[pixels[4 * (i * image_width + j) + k]];
                    for(int e = 0;e < 3;e++)
                    {
                        pixels[4*(i * image_width + j)+e] = result;
                    }
                }
            }
            return;
        }
        for(int i = 0 ;i < image_height;i ++)
        {
            for(int j = 0;j < image_width;j ++)
            {
                pixels[4 * (i * image_width + j) + k] = Histogram_result[pixels[4 * (i * image_width + j) + k]];
            }
        }
    }
}

bool InsideImage(int image_height, int image_width, int current_y, int current_x)
{
    return (current_y >= 0 && current_x >=0) && (current_y < image_height && current_x < image_width);
}

void SetSmoothingMatrix(int matrix_size, std::vector<std::vector<float>> &Matrix, float &matrix_scalar, float &matrix_scalar_DIV)
{
    if(matrix_size == 3)
    {
        matrix_scalar = 1;
        matrix_scalar_DIV = 10;
        for(int i = 0;i < matrix_size;i ++)
        {
            for(int j = 0;j < matrix_size;j ++)
            {
                Matrix[i][j] = 1;
            }
        }
        Matrix[1][1] = 2;
    }
    else if(matrix_size == 5)
    {
        matrix_scalar = 1;
        matrix_scalar_DIV = 273;
        Matrix[0][0] = 1;
        Matrix[0][1] = 4;
        Matrix[0][2] = 7;
        Matrix[0][3] = 4;
        Matrix[0][4] = 1;
        Matrix[1][0] = 4;
        Matrix[1][1] = 16;
        Matrix[1][2] = 26;
        Matrix[1][3] = 16;
        Matrix[1][4] = 4;
        Matrix[2][0] = 7;
        Matrix[2][1] = 26;
        Matrix[2][2] = 41;
        Matrix[2][3] = 26;
        Matrix[2][4] = 7;
        Matrix[3][0] = 4;
        Matrix[3][1] = 16;
        Matrix[3][2] = 26;
        Matrix[3][3] = 16;
        Matrix[3][4] = 4;
        Matrix[4][0] = 1;
        Matrix[4][1] = 4;
        Matrix[4][2] = 7;
        Matrix[4][3] = 4;
        Matrix[4][4] = 1;
    }
}

void SetEdgeDetectionMatrix(int matrix_size, std::vector<std::vector<float>> &Matrix, float &matrix_scalar, float &matrix_scalar_DIV)
{
    if(matrix_size == 3)
    {
        matrix_scalar = 1;
        matrix_scalar_DIV = 1;
        Matrix[0][0] = 2;
        Matrix[0][1] = -1;
        Matrix[0][2] = 2;
        Matrix[1][0] = -1;
        Matrix[1][1] = -4;
        Matrix[1][2] = -1;
        Matrix[2][0] = 2;
        Matrix[2][1] = -1;
        Matrix[2][2] = 2;
    }
    else if(matrix_size == 5)
    {
        matrix_scalar = 1;
        matrix_scalar_DIV = 1;
        Matrix[0][0] = 0;
        Matrix[0][1] = 0;
        Matrix[0][2] = -1;
        Matrix[0][3] = 0;
        Matrix[0][4] = 0;
        Matrix[1][0] = 0;
        Matrix[1][1] = -1;
        Matrix[1][2] = -2;
        Matrix[1][3] = -2;
        Matrix[1][4] = 0;
        Matrix[2][0] = -1;
        Matrix[2][1] = -2;
        Matrix[2][2] = 16;
        Matrix[2][3] = -2;
        Matrix[2][4] = -1;
        Matrix[3][0] = 0;
        Matrix[3][1] = -1;
        Matrix[3][2] = -2;
        Matrix[3][3] = -2;
        Matrix[3][4] = 0;
        Matrix[4][0] = 0;
        Matrix[4][1] = 0;
        Matrix[4][2] = -1;
        Matrix[4][3] = 0;
        Matrix[4][4] = 0;
    }
}

void Convolution(int &image_width, int &image_height, GLubyte* pixels, int matrix_size, float matrix_scalar, std::vector<std::vector<float>> Matrix)
{
    for(int k = 0;k < 3;k++) // RGB
    {
        std::vector<float> output;
        for(int i = 0;i < image_height;i ++)
        {
            for(int j = 0;j < image_width;j ++)
            {
                int current_i = i - matrix_size/2, current_j = j - matrix_size/2;
                float result = 0.f;
                for(int M_i = 0; M_i < matrix_size; M_i++)
                {
                    for(int M_j = 0;M_j < matrix_size;M_j++)
                    {
                        if(InsideImage(image_height,image_width,current_i+M_i,current_j+M_j))
                        {
                            result += pixels[4*((current_i+M_i)*image_width+(current_j+M_j))+k] * Matrix[matrix_size-M_i-1][matrix_size-M_j-1] * matrix_scalar;
                        }
                    }
                }
                if(result > 255) result = 255;
                if(result < 0) result = 0;
                //pixels[4*(i * image_width + j)+k] = result;
                output.emplace_back(result);
            }
        }
        //printf("%d",output.size());
        for(int i = 0;i < image_height;i ++)
        {
            for(int j = 0;j < image_width;j ++)
            {
                pixels[4*(i * image_width + j)+k] = output[i*image_width+j];
            }
        }
        output.clear();
    }
}

void CalculateBrightnessHistogram(int &image_width, int &image_height, GLuint &image_texture, int grayscale_type, float histogram[])
{
    GLubyte* pixels = Load_Texture(image_width, image_height, image_texture);
    // greyscale
    if(grayscale_type == 0) AverageGrayScale(image_width,image_height,pixels);
    else if(grayscale_type == 1) WeightedGrayScale(image_width,image_height,pixels);
    //Save_Texture(image_width,image_height,image_texture,pixels);
    for(int i = 0;i < 256;i ++) histogram[i] = 0;
    for(int i = 0;i < image_height;i ++)
    {
        for(int j = 0;j < image_width;j ++)
        {
            histogram[pixels[4*(i*image_width+j)]] ++;
        }
    }
    //for(int i = 0;i < 255;i ++) printf("%f\n",histogram[i]);
    free(pixels);
}

void DrawBrightnessHistogram(float histogram[])
{
    ImGui::Text("Brightness Histogram:");
    float max_value = 0;
    for(int i = 0;i < 255;i ++)
    {
        if(histogram[i] > max_value) max_value = histogram[i];
    }
    float graph_width = ImGui::GetWindowSize().x-20;
    ImGui::PlotHistogram("##BrightnessHistogram", histogram, 256, 0, NULL, 0.0f, max_value, ImVec2(graph_width, 300.0f));
    ImDrawList* draw_list = ImGui::GetWindowDrawList();
    for(int i = 0;i < 256;i ++)
    {
        ImVec2 pos = ImGui::GetCursorScreenPos();
        ImVec2 marker_min = ImVec2(pos.x + graph_width/256*i, pos.y);
        ImVec2 marker_max = ImVec2(pos.x + graph_width/256*(i+1), pos.y + 20);
        draw_list->AddRectFilled(marker_min, marker_max, IM_COL32(i, i, i, 255));
    }
    //ImGui::SetCursorPos(ImVec2(ImGui::GetCursorPos().x,ImGui::GetCursorPos().y+15));
    //ImGui::Text("0");
}

float CalculateNoiseHistogram(float StandardDeviation, int range_type)
{
    float Range = 0;
    if(range_type == 1) Range = StandardDeviation*5;
    else if(range_type == 0) Range = 1;
    float GaussianDistributionBinsize = Range*2.f/GaussianDistributionBinCount;

    for(int i = 0;i < GaussianDistributionBinCount;i ++)
    {
        GaussianDistribution[i] = 0; // clear out Gaussian array
    }

    for(int i = 0;i < z.size();i ++)
    {
        int z_bin = (int)((z[i]+Range)/GaussianDistributionBinsize);
        if(z[i] >= -1*Range && z[i] <= 1*Range) GaussianDistribution[z_bin] ++;
    }
    return Range;
}

void DrawNoiseHistogram(float histogram[], float range)
{
    ImGui::Text("Noise Histogram:");
    ImGui::SameLine(ImGui::GetWindowWidth()-160); ImGui::Text("[-%f, %f]",range,range);
    float max_value = 0;
    for(int i = 0;i < GaussianDistributionBinCount;i ++)
    {
        if(histogram[i] > max_value) max_value = histogram[i];
    }
    float graph_width = ImGui::GetWindowSize().x-20;
    ImGui::PlotHistogram("##NoiseHistogram", histogram, GaussianDistributionBinCount, 0, NULL, 0.0f, max_value, ImVec2(graph_width, 300.0f));
    //ImDrawList* draw_list = ImGui::GetWindowDrawList();
}

int main( void )
{
    srand( time(NULL) ); // set up random

    ShowWindow( GetConsoleWindow(), SW_HIDE );

    // - - - Initialize GL, UI - - -
    // Initialize GLFW
    if( !glfwInit() )
    {
        fprintf( stderr, "Failed to initialize GLFW\n" );
        getchar();
        return -1;
    }

    glfwWindowHint(GLFW_RESIZABLE, GLFW_FALSE);
    glfwWindowHint(GLFW_SAMPLES, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // To make MacOS happy; should not be needed
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    // Open a window and create its OpenGL context
    Window = glfwCreateWindow( Window_x, Window_y, "AIP 61047001s", NULL, NULL);
    if( Window == NULL ){
        fprintf( stderr, "Failed to open GLFW window. If you have an Intel GPU, they are not 3.3 compatible. Try the 2.1 version of the tutorials.\n" );
        getchar();
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(Window);

    // Initialize GLEW
    if (glewInit() != GLEW_OK) {
        fprintf(stderr, "Failed to initialize GLEW\n");
        getchar();
        glfwTerminate();
        return -1;
    }

    // Ensure we can capture the escape key being pressed below
    glfwSetInputMode(Window, GLFW_STICKY_KEYS, GL_TRUE);

    // Setup Dear ImGui context
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    io.IniFilename = NULL; // disable imgui.ini generate
    io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;           // imgui Docking Enable
    //io.ConfigFlags |= ImGuiConfigFlags_ViewportsEnable;       // Multiple Viewport Enable
    //io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;     // Enable Keyboard Controls
    //io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;      // Enable Gamepad Controls

    // Setup Dear ImGui style
    ImGui::StyleColorsDark();
    //ImGui::StyleColorsClassic();

    // Setup Platform/Renderer backends
    ImGui_ImplGlfw_InitForOpenGL(Window, true);
    ImGui_ImplOpenGL3_Init("#version 130");

    // Dark blue background
    glClearColor(0.0f, 0.0f, 0.4f, 0.0f);

    ImGuiWindowFlags docking_window_flags = 0;
    docking_window_flags |= ImGuiWindowFlags_NoTitleBar;
    docking_window_flags |= ImGuiWindowFlags_NoMove;
    docking_window_flags |= ImGuiWindowFlags_NoResize;
    docking_window_flags |= ImGuiWindowFlags_NoBringToFrontOnFocus;

    ImGuiWindowFlags working_window_flags = 0;
    //working_window_flags |= ImGuiWindowFlags_NoTitleBar;
    working_window_flags |= ImGuiWindowFlags_NoMove;

    int tmp = 0;

    int current_input = 1; // select analysis image 0 -> input 1 -> output
    const char* analysis_input[] = { "Input  Image", "Output Image" }; // combo for analysis
    int current_grayscale = 0; // select GrayScale 0 -> Average 1 -> Weighted
    int haar_current_grayscale = 0; // select GrayScale 0 -> Average 1 -> Weighted
    const char* analysis_grayscale[] = { "Average", "Weighted" }; // combo for analysis
    int current_analysis = 0; //
    const char* analysis_outputanalysis[] = { "Brightness Analysis" ,"Noise Analysis"}; // combo for analysis
    const char* analysis_inputanalysis[] = { "Brightness Analysis" }; // combo for analysis
    int current_range = 0; //
    const char* analysis_range[] = { "Static  Range" ,"Dynamic Range"}; // combo for Range
    bool if_analysis = false;
    const char* action_matrixsize[] = { "3x3", "5x5" }; // combo for convolution matrix size
    int current_matrixsize = 0; //
    int MatrixSize=3;
    const char* action_filter[] = { "Smoothing", "Edge Detect", "Custom" }; // combo for convolution matrix size
    int current_filter = 0; //

    // convolution matrix
    std::vector<std::vector<float>> Matrix3x3(3, std::vector<float>(3,1));
    float Matrix3x3Scalar = 1, Matrix3x3ScalarDiv = 1;
    std::vector<std::vector<float>> Matrix5x5(5, std::vector<float>(5,1));
    float Matrix5x5Scalar = 1, Matrix5x5ScalarDiv = 1;

    float Brightness_Histogram[256];

    float whitenoise_SD = 50.f;
    for(int i = 0;i < GaussianDistributionBinCount;i ++) GaussianDistribution[i] = 0;

    int haar_level = 1;
    float haar_scale = 25;

    int resize_height = 0,resize_width = 0;

    // - - - end initialize GL, UI - - -

    int input_image_width = 0;
    int input_image_height = 0;
    int output_image_width = 0;
    int output_image_height = 0;
    GLuint input_image_texture, output_image_texture;
    //std::string InputFileName = std::string("test.png");
    std::string InputFileName = std::string("default.jpg");
    std::string OutputFileName = std::string("default.jpg");
    ImVec4 tint_col = ImVec4(1.0f, 1.0f, 1.0f, 1.0f);   // No tint
    ImVec4 border_col = ImVec4(1.0f, 1.0f, 1.0f, 0.5f); // 50% opaque white
    // - - - loading initial image - - -

    LoadImage(InputFileName,input_image_width,input_image_height,output_image_width,output_image_height,input_image_texture,output_image_texture);

    resize_height = input_image_height;
    resize_width = input_image_width;

    // - - - end loading initial image - - -

    // - - - main draw loop - - -
    do{
        glClear( GL_COLOR_BUFFER_BIT );
        glfwPollEvents();

        // Start the Dear ImGui frame
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        // - - - UI layout - - -
        static ImGuiID dockspaceID = 0;
        {
            // set up a docking space
            const ImGuiViewport* Docking_Window = ImGui::GetMainViewport();
            ImGui::SetNextWindowPos(ImVec2(Docking_Window->WorkPos.x,Docking_Window->WorkPos.y),ImGuiCond_Once);
            ImGui::SetNextWindowSize(ImVec2(Window_x, Window_y), ImGuiCond_Once);
            ImGui::Begin("Docking", NULL ,docking_window_flags);
            dockspaceID = ImGui::GetID("HUB_DockSpace");
            ImGui::DockSpace(dockspaceID);
            //if (ImGui::DockBuilderGetNode(dockspaceID) == nullptr) printf("test\n");
            if(tmp == 0) // check if first run
            { // building dock layout
                tmp ++;
                ImGui::DockBuilderRemoveNode(dockspaceID); // Clear out existing layout
                //ImGui::DockBuilderAddNode(dockspaceID); //viewport->Size); // Add empty node
                ImGui::DockBuilderAddNode(dockspaceID, ImGuiDockNodeFlags_DockSpace); // Add empty node
                ImGui::DockBuilderSetNodeSize(dockspaceID, ImVec2(Window_x, Window_y));

                ImGuiID dock_main_id = dockspaceID; // This variable will track the document node, however we are not using it here as we aren't docking anything into it.
                ImGuiID dock_id_right;
                ImGuiID dock_id_left;
                ImGui::DockBuilderSplitNode(dock_main_id, ImGuiDir_Left, 0.49f, &dock_id_left, &dock_id_right);

                ImGui::DockBuilderDockWindow("  Input  ", dock_id_left);
                ImGui::DockBuilderDockWindow("Analysis", dock_id_left);
                ImGui::DockBuilderDockWindow(" Action ", dock_id_left);
                //ImGui::DockBuilderDockWindow("Docking", dock_main_id);
                ImGui::DockBuilderDockWindow(" Output ", dock_id_right);
                ImGui::DockBuilderFinish(dockspaceID);
            }
            ImGui::End();

            { // input window
                bool if_load_image = false;
                ImGui::Begin("  Input  ",NULL, working_window_flags);
                static char buf[256];
                strcpy(buf, InputFileName.c_str());
                if(ImGui::InputText("##inputPath", buf, IM_ARRAYSIZE(buf)))
                {
                    InputFileName = std::string(buf);
                    //if_load_image = true;
                }
                ImGui::SameLine();
                if (ImGui::Button(" Load "))
                {
                    if_load_image = true;
                }
                ImGui::SameLine();
                if (ImGui::Button("Browse"))
                {
                    if(BrowseFile(InputFileName))
                    {
                        if_load_image = true;
                    }
                }
                ImGui::Separator();

                ImVec2 window_size = ImGui::GetWindowSize();
                float imageScalar = window_size.x / fmax(input_image_width, input_image_height);
                ImGui::SetCursorPos(ImVec2((window_size.x - input_image_width * imageScalar) * 0.5,ImGui::GetCursorPos().y));
                ImGui::Image((void*)(intptr_t)input_image_texture, ImVec2(input_image_width * imageScalar, input_image_height * imageScalar));

                ImGui::End();

                if(if_load_image)
                {
                    LoadImage(InputFileName,input_image_width,input_image_height,output_image_width,output_image_height,input_image_texture,output_image_texture);
                    if_analysis = false;
                    resize_height = input_image_height;
                    resize_width = input_image_width;
                }
            }

            { // output window
                bool if_save_image = false;
                ImGui::Begin(" Output ",NULL, working_window_flags);
                static char buf[256];
                strcpy(buf, OutputFileName.c_str());
                if(ImGui::InputText("##outputPath", buf, IM_ARRAYSIZE(buf)))
                {
                    OutputFileName = std::string(buf);
                    //if_save_image = true;
                }
                ImGui::SameLine();
                if (ImGui::Button(" Save "))
                {
                    if_save_image = true;
                }
                ImGui::SameLine();
                if (ImGui::Button("Browse"))
                {
                    if(BrowseFileSave(OutputFileName, output_image_width, output_image_height, output_image_texture))
                    {
                        if_save_image = true;
                    }
                }
                ImGui::Separator();
                ImVec2 window_size = ImGui::GetWindowSize();
                float imageScalar = window_size.x / fmax(output_image_width, output_image_height);
                ImGui::SetCursorPos(ImVec2((window_size.x - input_image_width * imageScalar) * 0.5,ImGui::GetCursorPos().y));
                ImGui::Image((void*)(intptr_t)output_image_texture, ImVec2(output_image_width * imageScalar, output_image_height * imageScalar));
                ImGui::End();

                if(if_save_image)
                {
                    SaveImage(OutputFileName,output_image_width,output_image_height,output_image_texture);
                }
            }

            { // analysis window
                ImGui::Begin("Analysis",NULL, working_window_flags);
                // select input or output image
                ImGui::SetNextItemWidth(ImGui::GetWindowSize().x/2);
                if(ImGui::Combo(" Select Image", &current_input, analysis_input, IM_ARRAYSIZE(analysis_input)))
                {
                    if_analysis = false;
                    current_analysis = 0;
                }
                ImGui::Separator();

                ImGui::SetNextItemWidth(ImGui::GetWindowSize().x/4*3);
                if(current_input == 0)
                    if(ImGui::Combo("##analysis", &current_analysis, analysis_inputanalysis, IM_ARRAYSIZE(analysis_inputanalysis)))
                        if_analysis = false;
                if(current_input == 1)
                    if(ImGui::Combo("##analysis", &current_analysis, analysis_outputanalysis, IM_ARRAYSIZE(analysis_outputanalysis)))
                        if_analysis = false;
                ImGui::SameLine();
                if (ImGui::Button("Analysis"))
                {
                    if_analysis = true;
                }
                ImGui::Separator();
                if (ImGui::CollapsingHeader("Analysis Setting"))
                {
                    if(current_analysis == 0)
                    {
                        // grayscale method
                        ImGui::SetNextItemWidth(ImGui::GetWindowSize().x/2);
                        if(ImGui::Combo(" GrayScale Method", &current_grayscale, analysis_grayscale, IM_ARRAYSIZE(analysis_grayscale))) if_analysis = false;
                        ImGui::SameLine(); HelpMarker("Select the Desire GrayScale Method\nAverage:  (R + G + B) / 3\nWeighted: 0.299R + 0.587G + 0.114B");
                        ImGui::Separator();
                    }
                    else if(current_analysis == 1)
                    {
                        ImGui::SetNextItemWidth(ImGui::GetWindowSize().x/2);
                        if(ImGui::Combo(" Distribution Range", &current_range, analysis_range, IM_ARRAYSIZE(analysis_range))) if_analysis = false;
                        ImGui::SameLine(); HelpMarker("Select the Desire Distribution Range\nStatic  Range: [-1, 1]\nDynamic Range: Decide by Standard Deviation");
                        ImGui::Separator();
                    }
                }
                // histogram
                if(if_analysis)
                {
                    if(current_analysis == 0)
                    {
                        if(current_input == 0) CalculateBrightnessHistogram(input_image_width,input_image_height,input_image_texture,current_grayscale,Brightness_Histogram);
                        else if(current_input == 1) CalculateBrightnessHistogram(output_image_width,output_image_height,output_image_texture,current_grayscale,Brightness_Histogram);
                        ImGui::SetCursorPos(ImVec2(ImGui::GetCursorPos().x,ImGui::GetWindowHeight()/3-30));
                        DrawBrightnessHistogram(Brightness_Histogram);
                    }
                    else if(current_analysis == 1)
                    {
                        float range = CalculateNoiseHistogram(whitenoise_SD/255 ,current_range);
                        ImGui::SetCursorPos(ImVec2(ImGui::GetCursorPos().x,ImGui::GetWindowHeight()/3-30));
                        DrawNoiseHistogram(GaussianDistribution, range);
                    }
                }
                // reload button
                ImGui::End();
            }

            { // action window
                ImGui::Begin(" Action ",NULL, working_window_flags);
                ImGui::Text("Apply effect to the output image");
                ImGui::Separator();
                if (ImGui::CollapsingHeader("Reload"))
                {
                    ImGui::Separator();
                    ImGui::Text("\tReload From Input Image");
                    ImGui::SameLine(ImGui::GetWindowWidth() - 80);
                    if (ImGui::Button(" Reload ##Reload"))
                    {
                        GLubyte* pixels = Load_Texture(input_image_width,input_image_height,input_image_texture);
                        output_image_width = input_image_width;
                        output_image_height = input_image_height;
                        Save_Texture(output_image_width,output_image_height,output_image_texture,pixels);
                        free(pixels);
                        resize_height = input_image_height;
                        resize_width = input_image_width;
                        if(current_input == 1) if_analysis = false;
                    }
                }
                ImGui::Separator();
                if (ImGui::CollapsingHeader("Resize"))
                {
                    ImGui::Separator();
                    ImGui::Text("\tResize");
                    ImGui::SameLine(ImGui::GetWindowWidth() - 80);
                    if (ImGui::Button(" Resize ##resize"))
                    {
                        GLubyte* pixels = Load_Texture(output_image_width,output_image_height,output_image_texture);
                        GLubyte* resize_pixels = (GLubyte*)malloc(resize_width*resize_height*4*sizeof(GLubyte));
                        stbir_resize_uint8(pixels, output_image_width, output_image_height, 0, resize_pixels, resize_width, resize_height, 0, 4);
                        output_image_height = resize_height;
                        output_image_width = resize_width;
                        Save_Texture(output_image_width,output_image_height,output_image_texture,resize_pixels);
                        free(pixels);
                        free(resize_pixels);
                        if(current_input == 1) if_analysis = false;
                    }
                    ImGui::Text("\t\tHeight");
                    ImGui::SameLine(); ImGui::SetNextItemWidth(ImGui::GetWindowSize().x/4);
                    ImGui::DragInt("##height", &resize_height, 0.5f, 0, +INT_MAX, "%.3f", ImGuiSliderFlags_None);
                    ImGui::SameLine(ImGui::GetWindowWidth()/5*2+22);
                    ImGui::Text("\t\tWidth");
                    ImGui::SameLine(); ImGui::SetNextItemWidth(ImGui::GetWindowSize().x/4);
                    ImGui::DragInt("##width", &resize_width, 0.5f, 0, +INT_MAX, "%.3f", ImGuiSliderFlags_None);
                }
                ImGui::Separator();
                //ImGui::SetNextItemOpen(true, ImGuiCond_Appearing);
                if (ImGui::CollapsingHeader("GrayScale"))
                {
                    ImGui::Separator();
                    ImGui::Text("\tAverage  GrayScale");
                    ImGui::SameLine(); HelpMarker("Average:  (R + G + B) / 3");
                    ImGui::SameLine(ImGui::GetWindowWidth()-75);
                    if (ImGui::Button(" Apply ##Average"))
                    {
                        GLubyte* pixels = Load_Texture(output_image_width,output_image_height,output_image_texture);
                        AverageGrayScale(output_image_width,output_image_height,pixels);
                        Save_Texture(output_image_width,output_image_height,output_image_texture,pixels);
                        free(pixels);
                        if(current_input == 1) if_analysis = false;
                    }
                    ImGui::Separator();
                    ImGui::Text("\tWeighted GrayScale");
                    ImGui::SameLine(); HelpMarker("Weighted: 0.299R + 0.587G + 0.114B");
                    ImGui::SameLine(ImGui::GetWindowWidth()-75);
                    if (ImGui::Button(" Apply ##Weighted"))
                    {
                        GLubyte* pixels = Load_Texture(output_image_width,output_image_height,output_image_texture);
                        WeightedGrayScale(output_image_width,output_image_height,pixels);
                        Save_Texture(output_image_width,output_image_height,output_image_texture,pixels);
                        free(pixels);
                        if(current_input == 1) if_analysis = false;
                    }
                }
                ImGui::Separator();
                if (ImGui::CollapsingHeader("Noise"))
                {
                    ImGui::Separator();
                    ImGui::Text("\tWhite Noise");
                    ImGui::SameLine(ImGui::GetWindowWidth()/3+15);
                    //ImGui::Separator();
                    //ImGui::SetCursorPos(ImVec2(ImGui::GetCursorPos().x + 25,ImGui::GetCursorPos().y));
                    ImGui::SetNextItemWidth(ImGui::GetWindowSize().x/4); ImGui::DragFloat("Standard Deviation", &whitenoise_SD, 0.5f, 0, +FLT_MAX, "%.3f", ImGuiSliderFlags_None);
                    ImGui::Separator();
                    ImGui::Text("\t\tRGB  Noise");
                    ImGui::SameLine(); HelpMarker("Apply Seperate Noise Value to RGB");
                    ImGui::SameLine(ImGui::GetWindowWidth() - 75);
                    if (ImGui::Button(" Apply ##WhiteNoiseRGB")) {
                        GLubyte *pixels = Load_Texture(output_image_width, output_image_height, output_image_texture);
                        //WeightedGrayScale(output_image_width, output_image_height, pixels);
                        WhiteNoiseRGB(output_image_width,output_image_height,pixels,whitenoise_SD/255);
                        Save_Texture(output_image_width, output_image_height, output_image_texture, pixels);
                        free(pixels);
                        if (current_input == 1) if_analysis = false;
                    }
                    ImGui::Separator();
                    ImGui::Text("\t\tGray Noise");
                    ImGui::SameLine(); HelpMarker("Apply Single Noise Value to All RGB");
                    ImGui::SameLine(ImGui::GetWindowWidth() - 75);
                    if (ImGui::Button(" Apply ##WhiteNoiseGray"))
                    {
                        GLubyte *pixels = Load_Texture(output_image_width, output_image_height, output_image_texture);
                        //WeightedGrayScale(output_image_width, output_image_height, pixels);
                        WhiteNoise(output_image_width,output_image_height,pixels,whitenoise_SD/255);
                        Save_Texture(output_image_width, output_image_height, output_image_texture, pixels);
                        free(pixels);
                        if (current_input == 1) if_analysis = false;
                    }
                }
                ImGui::Separator();
                if(ImGui::CollapsingHeader("Wavelet Transform"))
                {
                    ImGui::Separator();
                    ImGui::Text("\tHaar Wavelet");
                    ImGui::SameLine(ImGui::GetWindowWidth() - 75);
                    if (ImGui::Button(" Apply ##HaarWavelet"))
                    {
                        GLubyte* pixels = Load_Texture(output_image_width,output_image_height,output_image_texture);
                        if(haar_current_grayscale == 0) AverageGrayScale(output_image_width,output_image_height,pixels);
                        else if(haar_current_grayscale == 1) WeightedGrayScale(output_image_width,output_image_height,pixels);
                        HaarWavelet(output_image_width,output_image_height,pixels,haar_level,haar_scale);
                        Save_Texture(output_image_width,output_image_height,output_image_texture,pixels);
                        free(pixels);
                        if(current_input == 1) if_analysis = false;
                    }
                    ImGui::Text("\t\tGreyScale Method");
                    ImGui::SameLine(); HelpMarker("Select the Desire GrayScale Method\nAverage:  (R + G + B) / 3\nWeighted: 0.299R + 0.587G + 0.114B");
                    ImGui::SameLine(ImGui::GetWindowWidth()/2+48);
                    // grayscale method
                    ImGui::SetNextItemWidth(ImGui::GetWindowSize().x/3);
                    if(ImGui::Combo("##GrayScale Method", &haar_current_grayscale, analysis_grayscale, IM_ARRAYSIZE(analysis_grayscale))) if_analysis = false;
                    ImGui::Text("\t\tLevel");
                    ImGui::SameLine(); ImGui::SetNextItemWidth(ImGui::GetWindowSize().x/4);
                    ImGui::DragInt("##level", &haar_level, 0.05f, 0, std::min(log2f(output_image_height),log2f(output_image_width)), "%.3f", ImGuiSliderFlags_None);
                    ImGui::SameLine(ImGui::GetWindowWidth()/5*2+22);
                    ImGui::Text("\t\tScale");
                    ImGui::SameLine(); ImGui::SetNextItemWidth(ImGui::GetWindowSize().x/4);
                    ImGui::DragFloat("##Scale", &haar_scale, 0.5f, 0, +FLT_MAX, "%.3f", ImGuiSliderFlags_None);
                    ImGui::Separator();
                }
                ImGui::Separator();
                if(ImGui::CollapsingHeader("Histogram Equalization"))
                {
                    ImGui::Separator();
                    ImGui::Text("\tGrey Equalization");
                    ImGui::SameLine(ImGui::GetWindowWidth() - 75);
                    if (ImGui::Button(" Apply ##GreyEqualization"))
                    {
                        GLubyte* pixels = Load_Texture(output_image_width,output_image_height,output_image_texture);
                        HistogramEqualization(output_image_width,output_image_height,pixels,0);
                        Save_Texture(output_image_width,output_image_height,output_image_texture,pixels);
                        free(pixels);
                        if(current_input == 1) if_analysis = false;
                    }
                    ImGui::Separator();
                    ImGui::Text("\tRGB  Equalization");
                    ImGui::SameLine(ImGui::GetWindowWidth() - 75);
                    if (ImGui::Button(" Apply ##RGBEqualization"))
                    {
                        GLubyte* pixels = Load_Texture(output_image_width,output_image_height,output_image_texture);
                        HistogramEqualization(output_image_width,output_image_height,pixels,1);
                        Save_Texture(output_image_width,output_image_height,output_image_texture,pixels);
                        free(pixels);
                        if(current_input == 1) if_analysis = false;
                    }
                }
                ImGui::Separator();
                if(ImGui::CollapsingHeader("Convolution"))
                {
                    ImGui::Separator();
                    ImGui::SetNextItemWidth(ImGui::GetWindowSize().x/2);
                    if(ImGui::Combo(" Matrix Size##Matrix size", &current_matrixsize, action_matrixsize, IM_ARRAYSIZE(action_matrixsize)))
                    {
                        if(current_matrixsize == 0) MatrixSize = 3;
                        if(current_matrixsize == 1) MatrixSize = 5;
                    }
                    ImGui::SetNextItemWidth(ImGui::GetWindowSize().x/2);
                    ImGui::Combo(" Set Matrix##Matrix", &current_filter, action_filter, IM_ARRAYSIZE(action_filter));
                    ImGui::SameLine(ImGui::GetWindowWidth() - 75);
                    if (ImGui::Button(" Apply ##Convolution"))
                    {
                        if(current_filter == 0)
                        {
                            if(current_matrixsize == 0) SetSmoothingMatrix(MatrixSize,Matrix3x3,Matrix3x3Scalar,Matrix3x3ScalarDiv);
                            if(current_matrixsize == 1) SetSmoothingMatrix(MatrixSize,Matrix5x5,Matrix5x5Scalar,Matrix5x5ScalarDiv);
                        }
                        if(current_filter == 1)
                        {
                            if(current_matrixsize == 0) SetEdgeDetectionMatrix(MatrixSize,Matrix3x3,Matrix3x3Scalar,Matrix3x3ScalarDiv);
                            if(current_matrixsize == 1) SetEdgeDetectionMatrix(MatrixSize,Matrix5x5,Matrix5x5Scalar,Matrix5x5ScalarDiv);
                        }

                        GLubyte* pixels = Load_Texture(output_image_width,output_image_height,output_image_texture);
                        if(current_matrixsize == 0) Convolution(output_image_width,output_image_height,pixels,MatrixSize,(Matrix3x3Scalar/Matrix3x3ScalarDiv), Matrix3x3);
                        if(current_matrixsize == 1) Convolution(output_image_width,output_image_height,pixels,MatrixSize,(Matrix5x5Scalar/Matrix5x5ScalarDiv), Matrix5x5);
                        Save_Texture(output_image_width,output_image_height,output_image_texture,pixels);
                        free(pixels);
                        if(current_input == 1) if_analysis = false;
                    }
                    if(current_filter == 2)
                    {
                        ImGui::Separator();
                        if(ImGui::CollapsingHeader("Custom Matrix"))
                        {
                            ImGui::Text("Matrix:");
                            if(current_matrixsize == 0) // 3x3
                            {
                                char MatrixLabel[5] = "##ij";
                                for(int i = 0;i < 3;i ++)
                                {
                                    for(int j = 0;j < 3;j ++)
                                    {
                                        MatrixLabel[2] = '0'+i;
                                        MatrixLabel[3] = '0'+j;
                                        if(j != 0) ImGui::SameLine();
                                        ImGui::SetNextItemWidth((ImGui::GetWindowSize().x-15*2-5*2)/3);
                                        ImGui::DragFloat(MatrixLabel, &Matrix3x3[i][j], 0.5f, -FLT_MAX, +FLT_MAX, "%.3f", ImGuiSliderFlags_None);
                                    }
                                }
                            }
                            else if(current_matrixsize == 1) // 5x5
                            {
                                char MatrixLabel[5] = "##ij";
                                for(int i = 0;i < 5;i ++)
                                {
                                    for(int j = 0;j < 5;j ++)
                                    {
                                        MatrixLabel[2] = '0'+i;
                                        MatrixLabel[3] = '0'+j;
                                        if(j != 0) ImGui::SameLine();
                                        ImGui::SetNextItemWidth((ImGui::GetWindowSize().x-15*2-5*5)/5);
                                        ImGui::DragFloat(MatrixLabel, &Matrix5x5[i][j], 0.5f, -FLT_MAX, +FLT_MAX, "%.3f", ImGuiSliderFlags_None);
                                    }
                                }
                            }
                            ImGui::Text("Scalar:"); ImGui::SameLine();
                            if(current_matrixsize == 0)
                            {
                                ImGui::SetNextItemWidth(ImGui::GetWindowSize().x/3.5);
                                ImGui::DragFloat("##MatrixScalar", &Matrix3x3Scalar, 0.5f, -FLT_MAX, +FLT_MAX, "%.3f", ImGuiSliderFlags_None);
                                ImGui::SameLine(); ImGui::Text("/");
                                ImGui::SameLine();
                                ImGui::SetNextItemWidth(ImGui::GetWindowSize().x/3.5);
                                ImGui::DragFloat("##MatrixScalarDiv", &Matrix3x3ScalarDiv, 0.5f, -FLT_MAX, +FLT_MAX, "%.3f", ImGuiSliderFlags_None);
                            }
                            else if(current_matrixsize == 1)
                            {
                                ImGui::SetNextItemWidth(ImGui::GetWindowSize().x/3.5);
                                ImGui::DragFloat("##MatrixScalar", &Matrix5x5Scalar, 0.5f, -FLT_MAX, +FLT_MAX, "%.3f", ImGuiSliderFlags_None);
                                ImGui::SameLine(); ImGui::Text("/");
                                ImGui::SameLine();
                                ImGui::SetNextItemWidth(ImGui::GetWindowSize().x/3.5);
                                ImGui::DragFloat("##MatrixScalarDiv", &Matrix5x5ScalarDiv, 0.5f, -FLT_MAX, +FLT_MAX, "%.3f", ImGuiSliderFlags_None);
                            }
                            ImGui::SameLine();
                            if(ImGui::Button("Set Sum"))
                            {
                                float sum = 0;
                                for(int i = 0;i < MatrixSize;i ++)
                                {
                                    for(int j = 0;j < MatrixSize;j ++)
                                    {
                                        if(current_matrixsize == 0) sum += Matrix3x3[i][j];
                                        else if(current_matrixsize == 1) sum += Matrix5x5[i][j];
                                    }
                                }
                                if(current_matrixsize == 0) Matrix3x3ScalarDiv = sum;
                                else if(current_matrixsize == 1) Matrix5x5ScalarDiv = sum;
                            }
                        }
                    }
                }
                ImGui::Separator();
                ImGui::End();
            }

            if(tmp == 1) // set input window to bo focused
            {
                ImGui::SetWindowFocus("  Input  ");
                tmp++;
            }

        }
        // - - - end UI layout - - -

        // Rendering
        ImGui::Render();

        //ImGui::UpdatePlatformWindows();
        //ImGui::RenderPlatformWindowsDefault();

        int display_w, display_h;
        glfwGetFramebufferSize(Window, &display_w, &display_h);
        glViewport(0, 0, display_w, display_h);
        glClear(GL_COLOR_BUFFER_BIT);
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());


        // Swap buffers
        glfwSwapBuffers(Window);

    } // Check if the ESC key was pressed or the window was closed
    while( glfwGetKey(Window, GLFW_KEY_ESCAPE ) != GLFW_PRESS &&
           glfwWindowShouldClose(Window) == 0 );
    // - - - end main draw loop - - -

    // - - - terminate - - -
    // Cleanup

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    glDeleteTextures(1, &input_image_texture);
    glDeleteTextures(1, &output_image_texture);

    glfwDestroyWindow(Window);
    glfwTerminate();

    return 0;
}