#include <iostream>
#include <OpenGL/gl.h>
#include <GLUT/glut.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <cfloat>

using namespace std;

enum status {
    light_source,
    no_intesection,
    intersection
};
struct point
{
    float x;
    float y;
    float z;
};

struct color
{
    float r;
    float g;
    float b;
};
struct sphere
{
    float radious;
    point spherePosition;
    color ks;
    color kd;
    color ka;
    float n;
};
struct source
{
    point sourcePosition;
    color is;
    color id;
    color ia;
};


int dimensions[2];
color background;
color global;


point viewer_v = {0.0, 0.0, 1.0};
float a = 1, b=0.1, c=0.01;
const int MAX_STEP = 2;
float viewport_size = 15.0;
point starting_directions = {0.0, 0.0, -1.0f};

vector<sphere> spheres;
vector<source> sources;

//paraser bierze sobie pierwsze słowo a potem kolejne parametry z niego
void fileParser(string path){
    ifstream infile(path);
    string line, type;
    while(getline(infile, line)){
        istringstream iss(line);
        iss>>type;
        if(!type.compare("dimensions")){
            iss>>dimensions[0]>>dimensions[1];
        }
        else if(!type.compare("background"))
            iss>>background.r>>background.g>>background.b;
        else if(!type.compare("global"))
            iss>>global.r>>global.g>>global.b;
        else if(!type.compare("sphere")) {
            sphere helpSphere;
            iss>>helpSphere.radious>>helpSphere.spherePosition.x>>helpSphere.spherePosition.y>>helpSphere.spherePosition.z;
            iss>>helpSphere.ks.r>>helpSphere.ks.g>>helpSphere.ks.b;
            iss>>helpSphere.kd.r>>helpSphere.kd.g>>helpSphere.kd.b;
            iss>>helpSphere.ka.r>>helpSphere.ka.g>>helpSphere.ka.b;
            iss>>helpSphere.n;
            spheres.push_back(helpSphere);
        }
        else if(!type.compare("source")){
            source helpSource;
            iss>>helpSource.sourcePosition.x>>helpSource.sourcePosition.y>>helpSource.sourcePosition.z;
            iss>>helpSource.is.r>>helpSource.is.g>>helpSource.is.b;
            iss>>helpSource.id.r>>helpSource.id.g>>helpSource.id.b;
            iss>>helpSource.ia.r>>helpSource.ia.g>>helpSource.ia.b;
            sources.push_back(helpSource);
        }
    }
}



// składowe intensywności świecenia źródła światła powodującego
/*************************************************************************************/
// Funkcja oblicza punkt przecięcia promienia i powierzchni sfery
// Argument p jest punktem początkowym promienia a d wektorem opisującym
// kierunek biegu promienia
// Funkcja zwraca 1 jeśli promień przecina sferę, 0 gdy nie przecina.
/*************************************************************************************/

point SingleInterces(point spherePosition, float radius, point p, point d, status *status1){
    float a, b, c, delta, r;
    a = d.x*d.x + d.y*d.y + d.z*d.z;  //tutaj jest opisana jakaś tam sfera
    b = 2*((p.x-spherePosition.x)*d.x + (p.y-spherePosition.y)*d.y + (p.z-spherePosition.z)*d.z);
    c = (float) (pow(p.x - spherePosition.x, 2) + pow(p.y - spherePosition.y, 2) + pow(p.z - spherePosition.z, 2) - pow(radius, 2));

    delta = b*b-4*a*c;
    point point1;
    if (delta>=0)                              // jest co najmniej jeden punkt przecięcia
    {
        r = (-b - sqrt(delta))/(2*a);     // parametr dla bliższego punktu przecięcia

        point1.x = p.x + r*d.x;    // współrzędne punktu przecięcia
        point1.y = p.y + r*d.y;
        point1.z = p.z + r*d.z;

        *status1 = intersection;                       // jest punkt przecięcia
    }
    else                                    // promień nie przecina sfery
        *status1 = no_intesection;
    return point1;
}


//ta funkcja normalizuje wektor
point Normalize(point p)
{
    float d =0.0;
    int i;
    d= (p.x*p.x)+(p.y*p.y)+(p.z*p.z);

    d=sqrt(d);

    if (d>0.0){
        p.x/=d;
        p.y/=d;
        p.z/=d;
    }
    return p;
}
// ta funkcja tworzy wektor normalny do powierzchnii danej sfery w punkcie q;
point Normal(point q, int sphereNumber){
    point returningPoint;
    returningPoint.x = q.x - spheres[sphereNumber].spherePosition.x;
    returningPoint.y = q.y - spheres[sphereNumber].spherePosition.y;
    returningPoint.z = q.z - spheres[sphereNumber].spherePosition.z;

    return Normalize(returningPoint);
}

/*************************************************************************************/
// Funkcja oblicza iloczyn skalarny wektorów
/*************************************************************************************/
float dotProduct(point p1, point p2)
{
    return (p1.x*p2.x+p1.y*p2.y+p1.z*p2.z);
}

//ta funkcja liczy sobie odpity promień - tak jak na stronie jarnickiego
point Reflect(point n, point light_vec){
    float n_dot_l = dotProduct(light_vec, n);
    point reflection_vector;

    reflection_vector.x = 2*(n_dot_l)*n.x - light_vec.x;
    reflection_vector.y = 2*(n_dot_l)*n.y - light_vec.y;
    reflection_vector.z = 2*(n_dot_l)*n.z - light_vec.z;

    return Normalize(reflection_vector);
}
float getLength(point vector){
    return sqrt(pow(vector.x,2)+pow(vector.y,2)+pow(vector.z,2));
}

//liczy odległość pomędzy wektorami
float calculateDistance(point source, point q){
    point vector;
    vector.x = source.x - q.x;
    vector.y = source.y - q.y;
    vector.z = source.z - q.z;
    return getLength(vector);

}

//tutaj sprawidzmy punkty przecięcia ze wszystkimi możliwymi punktami i znajdziemy najbliższy możliwy
point Intersect(point p, point d, status *status1, int *sphereNumber)
{
    float distance = FLT_MAX;
    point returningInter = {0, 0, 0};
    for(int i = 0; i<spheres.size(); ++i){
        status returningStatus;
        point sInter = SingleInterces(spheres[i].spherePosition, spheres[i].radious, p, d, &returningStatus);
        if(returningStatus&&calculateDistance(p, sInter)<distance) {
            distance = calculateDistance(p, sInter);
            *sphereNumber = i;
            *status1 = returningStatus;
            returningInter = sInter;
        }
    }

    return returningInter;
}

//model phonga dla pojedynczego źródła światła
color SinglePhong(int sphereNumber, int sourceNumber, point p, point n, point d){
    color returningColor;
    point light_vec;
    //  float d = calculateDistance(sources[sourceNumber].sourcePosition, q);
    light_vec.x = sources[sourceNumber].sourcePosition.x - p.x; // wektor wskazujący kierunek źródła
    light_vec.y = sources[sourceNumber].sourcePosition.y - p.y;
    light_vec.z = sources[sourceNumber].sourcePosition.z - p.z;

    float nl = dotProduct(n, Normalize(light_vec));
    float rv = dotProduct(Reflect(n, light_vec), viewer_v);
    float abdcd = 1;
    color ka = spheres[sphereNumber].ka;
    color ia = sources[sourceNumber].ia;
    color kd = spheres[sphereNumber].kd;
    color id = sources[sourceNumber].id;
    color ks = spheres[sphereNumber].ks;
    color is = sources[sourceNumber].is;
    float myN = spheres[sphereNumber].n;
    returningColor.r = (ka.r*ia.r)+(abdcd*((kd.r*id.r*nl)+(ks.r*is.r*pow(rv, myN))));
    returningColor.g = (ka.g*ia.g)+(abdcd*((kd.g*id.g*nl)+(ks.g*is.g*pow(rv, myN))));
    returningColor.b = (ka.b*ia.b)+(abdcd*((kd.b*id.b*nl)+(ks.b*is.b*pow(rv, myN))));
    returningColor.r/=2;
    returningColor.g/=2;
    returningColor.b/=2;


    return returningColor;
}

//ta funkcja wykonuje model phonga dla każdego możliwego źródła światła
color Phong(int sphereNumber, point light_vec, point n, point reflection, point q){
    color returningPhong = global;
    for(int i = 0; i<sources.size(); ++i){
        color phongColor = SinglePhong(sphereNumber, i ,q, n ,light_vec);
        returningPhong.r += phongColor.r*(1.f/(sources.size()+1));
        returningPhong.g += phongColor.g*(1.f/(sources.size()+1));
        returningPhong.b += phongColor.b*(1.f/(sources.size()+1));  //dodawanie kolorów odpywa sie na zasadzie liczenia średniej z kolorów dodawanych
    }

    return returningPhong;
}

//ta funkcja jest dokładnie opisana na stronie jarnickiego
color Trace(point p, point d, int *step){
    color local, reflected;
    point q;
    status myStatus;
    int sphereNumber;
    point n, r;
    if(*step > MAX_STEP)
        return background;

    q = Intersect(p, d, &myStatus, &sphereNumber);      //biere sobie sphere number żeby wiedzieć do która sfera jest najbliżej

    if(myStatus == no_intesection && *step == 0)
        return background;

    n = Normal(q, sphereNumber);
    r = Reflect(n, d);
   // local = Phong(sphereNumber, q, n, d);
    local = Phong(sphereNumber, d, n, r, q);
    *step+=1;
    reflected = Trace(q, r, step);

    color returningColor;
    returningColor.r = local.r+reflected.r;
    returningColor.g = local.g+reflected.g;
    returningColor.b = local.b+reflected.b;

    return returningColor;
}

void Display(void)
{
    float x_fl, y_fl;      // pozycja rysowanego piksela "zmiennoprzecinkowa"


    glClear(GL_COLOR_BUFFER_BIT);

    glFlush();

    // rysowanie pikseli od lewego górnego narożnika do prawego dolnego narożnika

    for (int y = dimensions[1]/2; y > -dimensions[1]/2; y--)
    {
        for (int x = -dimensions[0]/2; x < dimensions[0]/2; x++)
        {

            x_fl = (float)x/(dimensions[0]/viewport_size);
            y_fl = (float)y/(dimensions[1]/viewport_size);

            // przeliczenie pozycji(x,y) w pikselach na pozycję "zmiennoprzecinkową" w oknie obserwatora
            point starting_point;
            starting_point.x =  x_fl;
            starting_point.y =  y_fl;
            starting_point.z =  viewport_size;

            int step = 0;
            color pixel = Trace(starting_point,starting_directions, &step);
            if(step!=0) {
                pixel.r /= step;    //dzięki znanej liczbie kroków liczę sobie średnią z kolorów
                pixel.g /= step;
                pixel.b /= step;
            }
            glRasterPos3f(x_fl, y_fl, 0);

            // inkrementacja pozycji rastrowej dla rysowania piksela
            GLubyte pixel1[1][1][3];
            pixel1[0][0][0] = (GLubyte) (255 * pixel.r);
            pixel1[0][0][1] = (GLubyte) (255 * pixel.g);
            pixel1[0][0][2] = (GLubyte) (255 * pixel.b);
            glDrawPixels(1, 1, GL_RGB, GL_UNSIGNED_BYTE, pixel1);

            // narysowanie kolejnego piksela na ekranie
        }
    }
    glFlush();
}

/*************************************************************************************/
// Funkcja inicjalizująca definiująca sposób rzutowania
/*************************************************************************************/
void Myinit(void)
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-viewport_size/2, viewport_size/2, -viewport_size/2, viewport_size/2, -viewport_size/2, viewport_size/2);   //tak jak na stronie jarnickiego
    glMatrixMode(GL_MODELVIEW);
}

//to też ze strony
int main(int argc, char ** argv) {
    glutInit(&argc, argv);
    fileParser("scene.txt");
    cout<<global.r<<endl<<spheres[6].spherePosition.y;
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA);
    glutInitWindowSize(dimensions[0], dimensions[1]);
    glutCreateWindow("Ray Casting - oświetlona sfera");
    Myinit();
    glutDisplayFunc(Display);
    glutMainLoop();
    return 0;
}