#include  <iostream>
#include  <fstream>
#include  <string>
#include <vector>

using namespace std;
ifstream InFile;
ofstream OutFile;
double eps = 0.0000000000001;

class Vector3D {

public:
	double X;
	double Y; 
	double Z;
public:
	Vector3D(double dX, double dY, double dZ) {
		X = dX;
		Y = dY;
		Z = dZ;
	}

	Vector3D() {
		X = 0;
		Y = 0;
		Z = 0;
	}

	Vector3D(const Vector3D& A) {
		X = A.X;
		Y = A.Y;
		Z = A.Z;
	}

	Vector3D operator+(const Vector3D& vector2) const {
		Vector3D sum(this->X + vector2.X, this->Y + vector2.Y, this->Z + vector2.Z);
		return sum;
	}
	Vector3D operator-(const Vector3D& vector2) const {
		Vector3D sub(this->X - vector2.X, this->Y - vector2.Y, this->Z - vector2.Z);
		return sub;
	}

	Vector3D operator/(double x) const {
		Vector3D divide((this->X)/x, (this->Y)/x, (this->Z)/x);
		return divide;
	}

	Vector3D operator/(const Vector3D& vect) const {
		Vector3D divide((this->X) / vect.X, (this->Y) / vect.Y, (this->Z) / vect.Z);
		return divide;
	}

	bool operator==(const Vector3D& vector2) const{
		return (this->X == vector2.X && this->Y == vector2.Y && this->Z == vector2.Z);
	}

	
	bool check_zero() const{
		return (abs(this -> X) < eps || abs(this -> Y) < eps || abs(this->Z) < eps);
	}
};

class Segment3D {
public:
	Vector3D start;
	Vector3D end;

public:
	Segment3D(Vector3D vec1, Vector3D vec2) {
		start = vec1;
		end = vec2;
	}
};

bool is_skew(const Vector3D& line_1, const Vector3D& line_2, const Vector3D& line_3) {
	return(abs((line_1.X * line_2.Y * line_3.Z) + (line_2.X * line_3.Y * line_1.Z) + (line_1.Y * line_2.Z * line_3.X) -
		(line_1.Z * line_2.Y * line_3.X) - (line_1.X * line_3.Y * line_2.Z) - (line_2.X * line_1.Y * line_3.Z)) > eps);
}

bool is_colinear(const Vector3D& line_1, const Vector3D& line_2) {
	return(abs(line_1.Y * line_2.Z - line_1.Z * line_2.Y) < eps && abs(line_1.X * line_2.Y - line_1.Y * line_2.X) < eps &&
		abs(line_1.X * line_2.Z - line_1.Z * line_2.X) < eps);
}


//intersection of segments for identical lines

void intersect_identical(const Segment3D& seg_1, const Segment3D& seg_2, const Vector3D& line_1, const Vector3D& line_2) {
	if (abs(line_1.X) >= eps) {
		double lambda_2 = (seg_2.end.X - seg_1.start.X) / line_1.X;
		double lambda_1 = (seg_2.start.X - seg_1.start.X) / line_1.X;
		//one segment in another
		if (0 <= (lambda_1) && lambda_1 <=1 && 0 <= lambda_2 && lambda_2<=1) {
			OutFile << "Intersection: [(" << seg_2.start.X << ", " << seg_2.start.Y << "," << seg_2.start.Z << "), ("
				<< seg_2.end.X << ", " << seg_2.end.Y << ", " << seg_2.end.Z << ")]";
		}
		//end or start pf one segment in another
		else if (0 <= (lambda_1) && lambda_1 <= 1 && 1 <(lambda_2) ) {
			OutFile << "Intersection: [(" << seg_2.start.X << ", " << seg_2.start.Y << "," << seg_2.start.Z << "), ("
				<< seg_1.end.X << ", " << seg_1.end.Y << ", " << seg_1.end.Z << ")]";

		}
		else if (lambda_1 < 0 && 0 <= lambda_2 && lambda_2 <= 1) {
			OutFile << "Intersection: [(" << seg_1.start.X << ", " << seg_1.start.Y << "," << seg_1.start.Z << "), ("
				<< seg_2.end.X << ", " << seg_2.end.Y << ", " << seg_2.end.Z << ")]";
		}
		//one segment in another
		else if (lambda_1 < 0 && lambda_2>1) {
			OutFile << "Intersection: [(" << seg_1.start.X << ", " << seg_1.start.Y << "," << seg_1.start.Z << "), ("
				<< seg_1.end.X << ", " << seg_1.end.Y << ", " << seg_1.end.Z << ")]";
		}
		else if (lambda_1 > 1 || lambda_2 < 0) {
			OutFile << "No intersection";

		}
	}
	else if (abs(line_1.Y) >= eps) {
		double lambda_2 = (seg_2.end.Y - seg_1.start.Y) / line_1.Y;
		double lambda_1 = (seg_2.start.Y - seg_1.start.Y) / line_1.Y;
		if (0 <= (lambda_1) && lambda_1 <= 1 && 0 <= lambda_2 && lambda_2 <= 1) {
			OutFile << "Intersection: [(" << seg_2.start.X << ", " << seg_2.start.Y << "," << seg_2.start.Z << "), ("
				<< seg_2.end.X << ", " << seg_2.end.Y << ", " << seg_2.end.Z << ")]";
		}
		else if (0 <= (lambda_1) && lambda_1 <= 1 && 1 < (lambda_2)) {
			OutFile << "Intersection: [(" << seg_2.start.X << ", " << seg_2.start.Y << "," << seg_2.start.Z << "), ("
				<< seg_1.end.X << ", " << seg_1.end.Y << ", " << seg_1.end.Z << ")]";

		}
		else if (lambda_1 < 0 && 0 <= lambda_2 && lambda_2 <= 1) {
			OutFile << "Intersection: [(" << seg_1.start.X << ", " << seg_1.start.Y << "," << seg_1.start.Z << "), ("
				<< seg_2.end.X << ", " << seg_2.end.Y << ", " << seg_2.end.Z << ")]";
		}
		else if (lambda_1 < 0 && lambda_2>1) {
			OutFile << "Intersection: [(" << seg_1.start.X << ", " << seg_1.start.Y << "," << seg_1.start.Z << "), ("
				<< seg_1.end.X << ", " << seg_1.end.Y << ", " << seg_1.end.Z << ")]";
		}
		else if (lambda_1 > 1 || lambda_2 < 0) {
			OutFile << "No intersection";
		}
	}
	else if (abs(line_1.Z) >= eps) {
		double lambda_2 = (seg_2.end.Z - seg_1.start.Z) / line_1.Z;
		double lambda_1 = (seg_2.start.Z - seg_1.start.Z) / line_1.Z;
		if (0 <= (lambda_1) && lambda_1 <= 1 && 0 <= lambda_2 && lambda_2 <= 1) {
			OutFile << "Intersection: [(" << seg_2.start.X << ", " << seg_2.start.Y << "," << seg_2.start.Z << "), ("
				<< seg_2.end.X << ", " << seg_2.end.Y << ", " << seg_2.end.Z << ")]";
		}
		else if (0 <= (lambda_1) && lambda_1 <= 1 && 1 < (lambda_2)) {
			OutFile << "Intersection: [(" << seg_2.start.X << ", " << seg_2.start.Y << "," << seg_2.start.Z << "), ("
				<< seg_1.end.X << ", " << seg_1.end.Y << ", " << seg_1.end.Z << ")]";


		}
		else if (lambda_1 < 0 && 0 <= lambda_2 && lambda_2 <= 1) {
			OutFile << "Intersection: [(" << seg_1.start.X << ", " << seg_1.start.Y << "," << seg_1.start.Z << "), ("
				<< seg_2.end.X << ", " << seg_2.end.Y << ", " << seg_2.end.Z << ")]";
		}
		else if (lambda_1 < 0 && lambda_2>1) {
			OutFile << "Intersection: [(" << seg_1.start.X << ", " << seg_1.start.Y << "," << seg_1.start.Z << "), ("
				<< seg_1.end.X << ", " << seg_1.end.Y << ", " << seg_1.end.Z << ")]";
		}
		else if (lambda_1 > 1 || lambda_2 < 0) {
			OutFile << "No intersection";
		}
	}

}

void check_parallel_or_identical(const Segment3D& seg_1, const Segment3D& seg_2, const Vector3D& line_1, const Vector3D& line_2) {
	//if one of two vectors have all non-zero coordinates
	if (!line_1.check_zero() || !line_2.check_zero()) {
		if (!line_1.check_zero()) {
			Vector3D lambda = (seg_2.end - seg_1.start) / line_1;
			if (lambda.X == lambda.Y && lambda.Y == lambda.Z) {
				intersect_identical(seg_1, seg_2, line_1, line_2);
			}
			else {
				OutFile << "Segments are parallel and don't intersect\n";
			}
		}
		else if (!line_2.check_zero()) {
			Vector3D lambda = (seg_1.end - seg_2.start) / line_2;
			if (lambda.X == lambda.Y && lambda.Y == lambda.Z) {
				intersect_identical(seg_1, seg_2, line_1, line_2);
			}
			else {
				OutFile << "Segments are parallel and don't intersect\n";
			}

		}

	}
	else {
		//If two coordinates of vector  ~0
		if (abs(line_1.Z) < eps && abs(line_1.Y) < eps) {
			if (seg_1.start.Z == seg_2.start.Z && seg_1.start.Y == seg_2.start.Y) {
				intersect_identical(seg_1, seg_2, line_1, line_2);
			}
			else {
				OutFile << "Segments are parallel and don't intersect\n";
			}

		}
		else if (abs(line_1.X) < eps && abs(line_1.Y) < eps) {
			if (seg_1.start.X == seg_2.start.X && seg_1.start.Y == seg_2.start.Y) {
				intersect_identical(seg_1, seg_2, line_1, line_2);
			}
			else {
				OutFile << "Segments are parallel and don't intersect\n";
			}
		}
		else if (abs(line_1.Z) < eps && abs(line_1.X) < eps) {
			if (seg_1.start.X == seg_2.start.X && seg_1.start.Z == seg_2.start.Z) {
				intersect_identical(seg_1, seg_2, line_1, line_2);
			}
			else {
				OutFile << "Segments are parallel and don't intersect\n";
			}
		} //if one coordinate of vector ~0
		else if (abs(line_1.X) < eps) {
			if (seg_1.start.X == seg_2.start.X) {
				double l_y = (seg_2.end.Y - seg_1.start.Y) / (seg_1.end.Y - seg_1.start.Y);
				double l_z = (seg_2.end.Z - seg_1.start.Z) / (seg_1.end.Z - seg_1.start.Z);
				if (l_z == l_y) {
					intersect_identical(seg_1, seg_2, line_1, line_2);
				}
				else {
					OutFile << "Segments are parallel and don't intersect\n";
				}
			}
			else {
				OutFile << "Segments are parallel and don't intersect\n";
			}

		}
		else if (abs(line_1.Y) < eps) {
			if (seg_1.start.Y == seg_2.start.Y) {
				double l_x = (seg_2.end.X - seg_1.start.X) / (seg_1.end.X - seg_1.start.X);
				double l_z = (seg_2.end.Z - seg_1.start.Z) / (seg_1.end.Z - seg_1.start.Z);
				if (l_z == l_x) {
					intersect_identical(seg_1, seg_2, line_1, line_2);
				}
				else {
					OutFile << "Segments are parallel and don't intersect\n";
				}
			}
			else {
				OutFile << "Segments are parallel and don't intersect\n";
			}

		}
		else {
			if (seg_1.start.Z == seg_2.start.Z) {
				double l_x = (seg_2.end.X - seg_1.start.X) / (seg_1.end.X - seg_1.start.X);
				double l_y = (seg_2.end.Y - seg_1.start.Y) / (seg_1.end.Y - seg_1.start.Y);
				if (l_y == l_x) {
					intersect_identical(seg_1, seg_2, line_1, line_2);
				}
				else {
					OutFile << "Segments are parallel and don't intersect\n";
				}
			}
			else {
				OutFile << "Segments are parallel and don't intersect\n";
			}
		}
	}

}


void intersection_one_point(const Segment3D& seg_1, const Segment3D& seg_2, const Vector3D& line_1, const Vector3D& line_2) {
	vector<vector<double>> matrix(3, vector<double>(2));
	matrix[0][0] = line_1.X;
	matrix[0][1] = -line_2.X;
	matrix[1][0] = line_1.Y;
	matrix[1][1] = -line_2.Y;
	matrix[2][0] = line_1.Z;
	matrix[2][1] = -line_2.Z;
	double s = 1.5;
	double t = 1.5;
	vector<double> vector_right(3);
	vector_right[0] = seg_2.start.X - seg_1.start.X;
	vector_right[1] = seg_2.start.Y - seg_1.start.Y;
	vector_right[2] = seg_2.start.Z - seg_1.start.Z;
	if (abs(matrix[0][0])<eps) {
		if (abs(matrix[1][0]) > eps) {
			vector<double> change(2);
			for (int i = 0; i < 2; i++) {
				change[i] = matrix[0][i];
				matrix[0][i] = matrix[1][i];
				matrix[1][i] = change[i];
			}
			double change_vec = vector_right[0];
			vector_right[0] = vector_right[1];
			vector_right[1] = change_vec;



		}
		else {
			vector<double> change(2);
			for (int i = 0; i < 2; i++) {
				change[i] = matrix[0][i];
				matrix[0][i] = matrix[2][i];
				matrix[2][i] = change[i];

			}
			double change_vec = vector_right[0];
			vector_right[0] = vector_right[2];
			vector_right[2] = change_vec;


		}
	}

	for (int i = 1; i < 3; i++) {
		matrix[i][1] = matrix[i][1] - (matrix[i][0] * matrix[0][1]) / matrix[0][0];
		vector_right[i] = vector_right[i] - (matrix[i][0] * vector_right[0]) / matrix[0][0];
		matrix[i][0] = matrix[i][0] - matrix[i][0] * matrix[0][0] / matrix[0][0];
	}

	if (matrix[1][1] != 0) {
		s = vector_right[1] / matrix[1][1];
		t = (vector_right[0] - matrix[0][1] * s) / matrix[0][0];
		if (0 <= s && s <= 1 && 0 <= t && t <= 1) {
				OutFile << "X = " << t * line_1.X + seg_1.start.X << " Y = " << t * line_1.Y + seg_1.start.Y <<
					" Z = " << t * line_1.Z + seg_1.start.Z;
		}
		else {
				OutFile << "Segments don't intersect";
		}

	}
	else {
			s = vector_right[2] / matrix[2][1];
			t = (vector_right[0] - matrix[0][1] * s) / matrix[0][0];
			if (0 <= s && s <= 1 && 0 <= t && t <= 1) {
				OutFile << "X = " << t * line_1.X + seg_1.start.X << " Y = " << t * line_1.Y + seg_1.start.Y <<
					" Z = " << t * line_1.Z + seg_1.start.Z;
			}
			else {
				OutFile << "Segments don't intesect";
			}
	}


	

}

//function for intersection
void Intersect(Segment3D& seg_1, Segment3D& seg_2) {
	if (seg_1.start == seg_1.end || seg_2.start == seg_2.end) {
		OutFile << "Some segment is a point";

	}
	else if (seg_1.start == seg_2.start && seg_1.end == seg_2.end) {
		OutFile << "Segments are identical\n";
		OutFile << "Intesection [(" << seg_1.start.X << "," << seg_1.start.Y << "," << seg_1.start.Z <<
			"), (" << seg_1.end.X << "," << seg_1.end.Y << "," << seg_1.end.Z << ")]";
	}
	else {
		Vector3D line_1 = seg_1.end - seg_1.start;
		Vector3D line_2 = seg_2.end - seg_2.start;
		Vector3D line_3 = ((seg_2.end + seg_2.start) / 2) - ((seg_1.end + seg_1.start) / 2);
		//check tripple product
		if (is_skew(line_1, line_2, line_3)) {
			OutFile << "Segments are skrew and don't intersect" << "\n";
		}
		else {
			//check vector product 
			if (is_colinear(line_1, line_2)) {
				OutFile << "segments are collinear" << "\n";
				//check parallel or identical lines
				check_parallel_or_identical(seg_1, seg_2, line_1, line_2);
			}
			else {
				//Analysis of one point of intersection
				OutFile << "Segments may have one common point\n";
				if (seg_1.start == seg_2.end || seg_1.start == seg_2.start) {
					OutFile << "Intersection of segments (" << seg_1.start.X << "," << seg_1.start.Y <<
						"," << seg_1.start.Z << ")";

				}
				else if (seg_1.end == seg_2.end || seg_1.end == seg_2.start) {
					OutFile << "Intersection of segments (" << seg_1.end.X << "," << seg_1.end.Y <<
						"," << seg_1.end.Z << ")";

				}
				else {
					intersection_one_point(seg_1, seg_2, line_1, line_2);
				}
			}
		}
	}
	
}



int main() {
	InFile.open("input.txt");
	OutFile.open("output.txt");
	std::vector<double>points_1;
	std::vector<double>points_2;
	double point;
	for (int i = 0; i < 6; i++) {
		InFile >> point;
		points_1.push_back(point);
	}
	for (int i = 0; i < 6; i++) {
		InFile >> point;
		points_2.push_back(point);
	}
	InFile.close();

	Segment3D seg_1(Vector3D(points_1[0], points_1[1], points_1[2]), Vector3D(points_1[3], points_1[4], points_1[5]));
	Segment3D seg_2(Vector3D(points_2[0], points_2[1], points_2[2]), Vector3D(points_2[3], points_2[4], points_2[5]));
	Intersect(seg_1, seg_2);
	OutFile.close();
	return(0);
}