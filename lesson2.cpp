#include<cmath>
#include"tgaimage.h"
#include"geometry.h"
#include"model.h"
const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);
const TGAColor green = TGAColor(0,   255, 0,   255);
bool color_change = true;

int width = 1000;
int height = 1000;


double** allocateMatrix(int row, int col)
{
    double** matrix = new double* [row];
    for (int i = 0; i < row; ++i) {
        matrix[i] = new double[col] {0};
    }
    return matrix;
}

int deallocateMatrix(double** matrix, int row)
{
    for (int i = 0; i < row; ++i) {
        delete matrix[i];
    }
    delete[] matrix;
    return 0;
}
double** multiplyMatrix(double m1[3][3], int row1, int col1, double m2[3][1], int row2, int col2)
{
    if (col1 != row2)
        return nullptr;

    auto ret = allocateMatrix(row1, col2);

    int i, j, k;

    for (i = 0; i < row1; i++) {
        for (j = 0; j < col2; j++) {
            for (k = 0; k < col1; k++) {
                ret[i][j] += m1[i][k] * m2[k][j];
            }
        }
    }

    return ret;
}

void invdiagMat(int NbElement, double Mat[3][3])
{
    double** temp;
    int i;
    temp = new double* [NbElement];
    for (i = 0; i < NbElement; i++){
        temp[i] = new double[NbElement];
        for (int j = 0; j < NbElement; j++){
            temp[i][j] = 0;
        }
    }

    for (i = 0; i < NbElement; i++){
        for (int j = 0; j < NbElement; j++){
            temp[i][i] = 1 / Mat[i][i];
            if (j != i){
                temp[i][j] = -Mat[i][j] / Mat[i][i];
            }
            for (int k = 0; k < NbElement; k++){

                if (k != i){
                    temp[k][i] = Mat[k][i] / Mat[i][i];
                }
                if (j != i && k != i){
                    temp[k][j] = Mat[k][j] - Mat[i][j] * Mat[k][i] / Mat[i][i];
                }
            }
        }
        for (int i = 0; i < NbElement; i++){
            for (int j = 0; j < NbElement; j++){
                Mat[i][j] = temp[i][j];
            }
        }
    }
}

void line(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color) { 
    bool steep = false; 
    if (std::abs(x0-x1)<std::abs(y0-y1)) { // if the line is steep, we transpose the image 
        std::swap(x0, y0); 
        std::swap(x1, y1); 
        steep = true; 
    } 
    if (x0>x1) { // make it left−to−right 
        std::swap(x0, x1); 
        std::swap(y0, y1); 
    } 
    for (int x=x0; x<=x1; x++) { 
        float t = (x-x0)/(float)(x1-x0); 
        int y = y0*(1.-t) + y1*t; 
        if (steep) { 
            image.set(y, x, color); // if transposed, de−transpose 
        } else { 
            image.set(x, y, color); 
        } 
    } 
}


void line(Vec2i t0, Vec2i t1, TGAImage &image, TGAColor color){
    line(t0.x, t0.y, t1.x,t1.y,  image, color);
}

void triangle(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage &image, TGAColor color) {
    line(t0, t1, image, color); 
    line(t1, t2, image, color); 
    line(t2, t0, image, color); 
}

bool equal3Vec(Vec2i t0, Vec2i t1, Vec2i t2) {
    if(t0.x == t1.x && t1.x == t2.x && t0.y == t1.y && t1.y == t2.y)
        return true;
    return false;
}
void recentre(Vec2i t0, Vec2i t1, Vec2i t2,Vec2i barr,int tab[6]) {


    if (t0.x> barr.x && color_change) tab[0] = t0.x - 1;
    else if (t0.x< barr.x && color_change) tab[0] = t0.x + 1;
    else tab[0] = t0.x;

    if (t0.y> barr.y && !color_change) tab[1] = t0.y - 1;
        else if (t0.y< barr.y && !color_change) tab[1] = t0.y + 1;
        else tab[1] = t0.y;



    if (t1.x> barr.x && !color_change) tab[2] = t1.x - 1;
        else if (t1.x< barr.x && !color_change) tab[2] = t1.x + 1;
        else tab[2] = t1.x;

    if (t1.y> barr.y && color_change) tab[3] = t1.y - 1;
        else if (t1.y< barr.y && color_change) tab[3] = t1.y + 1;
        else tab[3] = t1.y;



    if (t2.x> barr.x && color_change) tab[4] = t2.x - 1;
            else if (t2.x< barr.x && color_change) tab[4] = t2.x + 1;
            else
                tab[4] = t2.x;

        if (t2.y> barr.y && !color_change) tab[5] = t2.y - 1;
            else if (t2.y< barr.y && !color_change) tab[5] = t2.y + 1;
            else
                tab[5] = t2.y;

}


void fill_trianglerec(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage& image, TGAColor color, Vec2i barr) {
    if (equal3Vec(t0, t1, t2))
        return;
    line(t0, t1, image, color); 
    line(t1, t2, image, color); 
    line(t2, t0, image, color);
    int tmp_value[6];
    recentre(t0,t1,t2,barr,tmp_value);
    Vec2i nt0(tmp_value[0], tmp_value[1]);
    triangle(nt0, t1, t2, image, color);
    Vec2i nt1(tmp_value[2], tmp_value[3]);
    triangle(t0, nt1, t2, image, color);
    Vec2i nt2(tmp_value[4], tmp_value[5]);
    triangle(t0, t1, nt2, image, color);
    if (color_change)
        color_change = false;
    else
        color_change = true;
    fill_trianglerec(nt0, nt1, nt2, image, color, barr);

}/*
void fill_triangle(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage& image, TGAColor color) {
    Vec2i Barricentre((t0.x + t1.x + t2.x) / 3,(t0.y + t1.y + t2.y) / 3 );
    
    fill_trianglerec(t0, t1, t2, image, color, Barricentre);
}*/

void fill_triangle(Vec3i t0, Vec3i t1, Vec3i t2, TGAImage& image, TGAColor color, float *zbuffer) {

    double matrice[3][3] = { {t0.x,t1.x,t2.x},{t0.y,t1.y,t2.y},{1,1,1} };

    invdiagMat(3,matrice);/*
    for (int i : {0, 1, 2}) {
        for (int j : {0, 1, 2})
            std::cerr << matrice[i][j] << " ";
        std::cerr << std::endl;
    }*/
    double vector[3][1] = { };
    for(int i = std::min(t0.x,std::min(t1.x,t2.x));i< std::max(t0.x,std::max(t1.x,t2.x)); i++)
        for(int j = std::min(t0.y,std::min(t1.y,t2.y));j< std::max(t0.y,std::max(t1.y,t2.y)); j++){
            vector[0][0] = i;
            vector[1][0] = j;
            vector[2][0] = 1;
            double ** verif = multiplyMatrix(matrice, 3, 3, vector, 3, 1);
            double z = t0.z * verif[0][0] + t1.z * verif[1][0] + t2.z * verif[2][0];

            if (verif[0][0] > -0.001 && verif[1][0] > -0.001 && verif[2][0] > -0.001 && zbuffer[i + j * width]<z){
                image.set(i, j, color);
                
                zbuffer[i + j * width] = z;
            }
            deallocateMatrix(verif, 3);
        }
        
}

int main(void)
{

    float* zbuffer = new float[width * height];
    for (int i = width * height; i--; zbuffer[i] = -std::numeric_limits<float>::max());

    TGAImage image(width,height, TGAImage::RGB);
    Model *model = new Model("african_head.obj");
    printf("le nombre de face : %d\n",model->nfaces());
    for (int i = 0; i < model->nfaces(); i++) {
        std::vector<int> face = model->face(i);
        Vec3i screen_coords[3];
        for (int j = 0; j < 3; j++) {
            Vec3f world_coords = model->vert(face[j]);
            screen_coords[j] = Vec3i((world_coords.x + 1.) * width / 2., (world_coords.y + 1.) * height / 2., world_coords.z*1000);
        }
        Vec3f n = ((model->vert(face[1]) - model->vert(face[0])) ^ (model->vert(face[2]) - model->vert(face[0]))).normalize();
        Vec3f l = Vec3f(13, -22, 3).normalize();
        //Vec3f l = Vec3f(-1, 1, -5).normalize();

        int intensity = std::max(0.f, n * l)*255.f;
        int color_intensity1 = ((model->vert(face[0]).z - 1) / -2)*255;
        int color_intensity2 = ((model->vert(face[1]).z - 1) / -2)*255;
        int color_intensity3 = ((model->vert(face[2]).z - 1) / -2)*255;
        int color_intensity = std::min(color_intensity1, std::min(color_intensity2, color_intensity3));
        intensity *= 1.8;
        fill_triangle(screen_coords[0], screen_coords[1], screen_coords[2],  image, TGAColor(std::min( intensity*211/255,255) ,std::min(intensity*195/255,255),std::min(intensity*123/255,255) ),zbuffer);
    }

    image.write_tga_file("image.tga");
    return 0;
}