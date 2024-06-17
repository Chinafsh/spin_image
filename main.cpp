#include<iostream>
#include<string>
#include <cassert>  // for assert function
#include<fstream>
#include <sstream>
#include<iomanip>
#include<vector>
#include<math.h>
#include <cfloat> // for max, min of float 
#include<algorithm>
//gargoyle100K, bunny
const std::string filename = "model/bunny.obj";
float mesh_resolution;    // the median of all edge lengths
int rows;
int cols;
float alpha_maxf = FLT_MIN;
float beta_maxf = FLT_MIN;


struct Vector2f{
    float X, Y;

    Vector2f(float _X, float _Y){
        X = _X;
        Y = _Y;
    }
};

struct Vector3f{
    float X, Y, Z;

    // Default Constructor
	Vector3f()
	{
    	X = 0.0f;
		Y = 0.0f;
		Z = 0.0f;
	}
		// Variable Set Constructor
	Vector3f(float X_, float Y_, float Z_)
	{
	    X = X_;
		Y = Y_;
		Z = Z_;
	}
    Vector3f& operator= (const Vector3f& right)
    {
        if(this != &right)
        {
            X = right.X;
            Y = right.Y;
            Z = right.Z;
        }
        return *this;
    }
	// Addition Operator Overload
	Vector3f operator+(const Vector3f& right) const
	{
        Vector3f res;
        res.X = X + right.X;
        res.Y = Y + right.Y;
        res.Z = Z + right.Z;
		return res;
	}
    // Subtraction Operator Overload
	Vector3f operator-(const Vector3f& right) const
	{
        Vector3f res;
        res.X = X - right.X;
        res.Y = Y - right.Y;
        res.Z = Z - right.Z;
		return res;
	}
	// Float Multiplication Operator Overload
	Vector3f operator*(const float& other) const
	{
		return Vector3f(this->X * other, this->Y * other, this->Z * other);
	}
	// Float Division Operator Overload
	Vector3f operator/(const float& other) const
	{
		return Vector3f(this->X / other, this->Y / other, this->Z / other);
	}        
    Vector3f operator+=(const Vector3f& right) const
    {
        return Vector3f(this->X+right.X, this->Y + right.Y, this->Z + right.Z);
    }

};


Vector3f CrossV3f(const Vector3f a, const Vector3f b)
{
	return Vector3f(a.Y * b.Z - a.Z * b.Y,
		a.Z * b.X - a.X * b.Z,
		a.X * b.Y - a.Y * b.X);
}

float MagnitudeV3f(const Vector3f in)
{
	return (sqrtf(powf(in.X, 2) + powf(in.Y, 2) + powf(in.Z, 2)));
}
float Norm2(const Vector3f in)
{
	return powf(in.X, 2) + powf(in.Y, 2) + powf(in.Z, 2);
}
float DotV3(const Vector3f a, const Vector3f b)
{
	return (a.X * b.X) + (a.Y * b.Y) + (a.Z * b.Z);
}

float AngleBetweenV3(const Vector3f a, const Vector3f b)
{
	float angle = DotV3(a, b);
	angle /= (MagnitudeV3f(a) * MagnitudeV3f(b));
	return angle = acosf(angle);
}

Vector3f Normalized(Vector3f a)
{
    return a/(MagnitudeV3f(a)*1.0);
}


struct Oriented_point{
    Vector3f Position;
    Vector3f Normal;
    std::vector<std::vector<float>> spin_image;


    Oriented_point(){
        Position = Vector3f(0,0,0);
        Normal = Vector3f(0,0,0);
    }

};



void loadObjFile(const std::string& filename, std::vector<Vector3f>& vertices, std::vector<std::vector<int>>& faces) {
    std::ifstream file(filename);
    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string token;
        iss >> token;
        if (token == "v") {
            Vector3f vertex;
            iss >> vertex.X >> vertex.Y >> vertex.Z;

            vertices.push_back(vertex);
        } else if (token == "f") {
            std::vector<int> face;
            int vIndex;
            while (iss >> vIndex) {
                face.push_back(vIndex-1); // Obj file indexing starts from 1
            }
            face.resize(3);
            faces.push_back(face);
        }
    }
}

float computeMeshResolution(std::vector<Vector3f>& verts, std::vector<std::vector<int>>& faces)
{
    // Mesh resolution is defined as the average of all edge lengths in a mesh. (Or median..)
    // std::vector<float> lengths;
    float sum = 0.0f;
    for(auto face : faces)
    {
        Vector3f p0 = verts[face[0]];
        Vector3f p1 = verts[face[1]];
        Vector3f p2 = verts[face[2]];

        sum += MagnitudeV3f(p0-p1);
        sum += MagnitudeV3f(p0-p2);
        sum += MagnitudeV3f(p2-p1);

    }
    
    return sum / (2.0*3.0*faces.size());
    
    
}


void computeFacesNormal(const std::vector<Vector3f>& points ,const std::vector<std::vector<int>>& faces, std::vector<Vector3f>& Normals)
{
    for(auto face : faces){
        Vector3f a = points[face[0]];
        Vector3f b = points[face[1]];
        Vector3f c = points[face[2]];

        Vector3f a2b = b-a;
        Vector3f a2c = c-a;

        Vector3f normal = CrossV3f(a2b, a2c);
        normal = Normalized(normal);
        Normals.push_back(normal);
    }
}


void getOriented_points(std::vector<Vector3f>& verts,  std::vector<std::vector<int>>& faces, std::vector<Oriented_point>& Os)
{

    std::vector<Vector3f> normals;
    computeFacesNormal(verts,faces,normals);

    for(int i=0; i<verts.size(); i++){
        Oriented_point O;
        O.Position = Vector3f(0,0,0);
        O.Normal = Vector3f(0,0,0);
        Os.push_back(O);
    }
    Os.resize(verts.size());    
    for(int i=0; i<faces.size(); i++)
    {

        int p0 = faces[i][0];
        int p1 = faces[i][1];
        int p2 = faces[i][2];

        if(p0<0||p0>=Os.size()||p1<0||p1>=Os.size()||p2<0||p2>=Os.size()){
            std::cout << "failure!!" << std::endl;
        }
        Os[p0].Normal = Os[p0].Normal + normals[i];
        Os[p1].Normal = Os[p1].Normal + normals[i];

        Os[p2].Normal = Os[p2].Normal + normals[i];
        

    }

    for(int i=0; i<Os.size(); i++){
        Os[i].Position = verts[i];
        Os[i].Normal = Normalized(Os[i].Normal);
    }

    // save data points and its normal vector
    std::ofstream outFile("data.txt");
    if (outFile.is_open()) {
        for (size_t i = 0; i < Os.size(); ++i) {
            outFile << std::fixed << std::setprecision(6) <<
            Os[i].Position.X << " " << Os[i].Position.Y << " " << Os[i].Position.Z<< " "
                    << Os[i].Normal.X << " " << Os[i].Normal.Y<< " " << Os[i].Normal.Z << std::endl;
        }
        outFile.close();
        std::cout << "Data saved to data.txt" << std::endl;
    } else {
        std::cerr << "Unable to open file for writing" << std::endl;
    }

}


Vector2f SpinMap(Vector3f& Position, Oriented_point& O)
{
    float beta = DotV3(O.Normal, Position-O.Position);// signed
    float alpha = sqrtf(Norm2(Position-O.Position) - beta*beta);// positive
    return Vector2f(alpha, beta);
}   

void ComputeMaxAB(std::vector<Oriented_point>& Os, float& a_max, float& b_max)
{
    int n = Os.size();
    for(int i=0; i<n; i++){
        for(int j=i+1; j<n; j++){
            Vector2f ab = SpinMap(Os[j].Position, Os[i]);
            a_max = std::max(a_max, ab.X);
            b_max = std::max(b_max, fabs(ab.Y));
        }
    }
}

std::vector<int> SpinImageBin(Vector2f& ab)
{

    int i = std::floor((beta_maxf-ab.Y)/mesh_resolution);       // row -- beta
    int j = std::floor(ab.X/mesh_resolution);                   // col -- alpha
    std::vector<int> res = {i, j};
    return res;
}

Vector2f BilinearWeights(Vector2f& ab, std::vector<int>& ij)
{
    float a = ab.X - ij[0]*mesh_resolution;
    float b = ab.Y - ij[1]*mesh_resolution;
    return Vector2f(a,b);
}



void CreateSpinImages(Oriented_point& O, std::vector<Oriented_point>& Os)
{
    int cnt = 0;
    int rows = O.spin_image.size();
    int cols = O.spin_image[0].size();
    for(int i=0; i<Os.size(); i++){
        float angle = AngleBetweenV3(O.Normal, Os[i].Normal);
        if(angle > 0&&angle<M_PI/3){
            // cnt++;
            Vector2f alphabeta = SpinMap(Os[i].Position, O);

            std::vector<int> ij = SpinImageBin(alphabeta);

            Vector2f ab = BilinearWeights(alphabeta, ij);

            if(ij[0]>=0 && ij[0]<rows && ij[1]>=0 && ij[1] < cols)
                O.spin_image[ij[0]][ij[1]] += (1-ab.X)*(1-ab.Y);

            if(ij[0]+1>=0&&ij[0]+1<rows&&ij[1]>=0 && ij[1] < cols)
                O.spin_image[ij[0]+1][ij[1]] += ab.X*(1-ab.Y);

            if(ij[0]>=0 && ij[0]<rows && ij[1]+1>=0 && ij[1]+1 < cols)
                O.spin_image[ij[0]][ij[1]+1] += (1-ab.X)*ab.Y;

            if(ij[0]+1>=0 && ij[0]+1<rows && ij[1]+1>=0 && ij[1]+1 < cols)
                O.spin_image[ij[0]+1][ij[1]+1] += ab.X*ab.Y;
        }

    }
    // std::cout << "contributions  = "<< cnt << std::endl;
}





int main()
{




    std::vector<Vector3f> verts;
    std::vector<std::vector<int>> faces;
    
    loadObjFile(filename, verts, faces);
    //bin size = mesh resolution
    mesh_resolution = computeMeshResolution(verts,faces);
    // Debug
    std::cout << "vertex: " << verts.size() << std::endl;
    std::cout << "face: " << faces.size() << std::endl;
    std::cout << "mesh resolution: " << mesh_resolution << std::endl;

    std::vector<Oriented_point> Oriented_points(verts.size());
    
    getOriented_points(verts, faces, Oriented_points);
    
    // maximum size of the object in spin-map coordinates.
    ComputeMaxAB(Oriented_points, alpha_maxf, beta_maxf);
    

    // the size of the spin-image (rows, cols)
    int rows = 2*beta_maxf/mesh_resolution + 1;
    int cols = alpha_maxf/mesh_resolution + 1;
    std::cout << "Spin-Image size: (" << rows << ", " << cols << ")" << std::endl;

    for(int i=0; i<Oriented_points.size(); i++){
        Oriented_points[i].spin_image.resize(rows, std::vector<float>(cols, 0.0f));
    }
    
    for(int k=0; k<Oriented_points.size(); k++){
    CreateSpinImages(Oriented_points[k],Oriented_points);
        std::string name = "SpinImages/spinImage" + std::to_string(k) + ".txt";
        std::ofstream outFile(name);
        if (outFile.is_open()) {
            for (int i = 0; i < rows; i++) {
                for(int j = 0; j <cols; j++){
                    if(j==cols-1)
                        outFile << std::fixed << std::setprecision(6) << Oriented_points[k].spin_image[i][j] << std::endl;
                    else
                        outFile << std::fixed << std::setprecision(6) << Oriented_points[k].spin_image[i][j] <<" ";
                }
            }
            outFile.close();
            // std::cout << "Data saved to spinImage.txt" << std::endl;
        } else {
            std::cerr << "Unable to open file for writing" << std::endl;
        }
    }
    std::cout << "Data saved to spinImages file" << std::endl;

    return 0;
}


