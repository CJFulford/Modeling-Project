#include "Header.h"


void readFromObj(
	const string fileDir,
	vector<MyLight>&  lightVec,
	vector<MySphere>& sphereVec,
	vector<MyPlane>& planeVec,
	vector<MyTriangle>& triangleVec
)
{
	ifstream file(fileDir);
	if (file.is_open())
	{
		string line;
		while (getline(file, line))
		{
			switch (line[0])
			{
				case 'v':
					switch (line[1])
					{
					case ' ':
						break;
					case 't':
						break;
					case 'n':
						break;
					case 'p':
						break;
					default:
						break;
					}
					break;
				case 'f':
					break;
				default:
					continue;
					break;
			}

			if (line.compare(0, string::npos, "light {") == 0)
			{
				MyLight light;
				getline(file, line);
				getPointFromLine(line, &light.point);

				getline(file, line);
				getPointFromLine(line, &light.colour);

				light.ambient = light.colour * 0.6f;

				lightVec.push_back(light);

			}

			else if (line.compare(0, string::npos, "sphere {") == 0)
			{
				MySphere sphere;
				getline(file, line);
				getPointFromLine(line, &sphere.center);

				getline(file, line);
				sphere.radius = stof(line.substr(line.find_first_not_of(' '), string::npos));

				getline(file, line);
				getPointFromLine(line, &sphere.colour);

				getline(file, line);
				getPointFromLine(line, &sphere.specular);

				getline(file, line);
				sphere.phong = stof(line.substr(line.find_first_not_of(' '), string::npos));

				getline(file, line);
				sphere.reflect = stof(line.substr(line.find_first_not_of(' '), string::npos));

				sphereVec.push_back(sphere);
			}

			else if (line.compare(0, string::npos, "plane {") == 0)
			{
				MyPlane plane;
				getline(file, line);
				getPointFromLine(line, &plane.normal);

				getline(file, line);
				getPointFromLine(line, &plane.point);

				getline(file, line);
				getPointFromLine(line, &plane.colour);

				getline(file, line);
				getPointFromLine(line, &plane.specular);

				getline(file, line);
				plane.phong = stof(line.substr(line.find_first_not_of(' '), string::npos));

				getline(file, line);
				plane.reflect = stof(line.substr(line.find_first_not_of(' '), string::npos));

				planeVec.push_back(plane);
			}

			else if (line.compare(0, string::npos, "triangle {") == 0)
			{
				MyTriangle tri;
				getline(file, line);
				getPointFromLine(line, &tri.p1);

				getline(file, line);
				getPointFromLine(line, &tri.p2);

				getline(file, line);
				getPointFromLine(line, &tri.p3);

				getline(file, line);
				getPointFromLine(line, &tri.colour);

				getline(file, line);
				getPointFromLine(line, &tri.specular);

				getline(file, line);
				tri.phong = stof(line.substr(line.find_first_not_of(' '), string::npos));

				getline(file, line);
				tri.reflect = stof(line.substr(line.find_first_not_of(' '), string::npos));

				tri.a = tri.p1.x - tri.p2.x;
				tri.b = tri.p1.y - tri.p2.y;
				tri.c = tri.p1.z - tri.p2.z;

				tri.d = tri.p1.x - tri.p3.x;
				tri.e = tri.p1.y - tri.p3.y;
				tri.f = tri.p1.z - tri.p3.z;

				triangleVec.push_back(tri);
			}
		}
	}
	else cout << "didnt open" << endl;
	file.close();
}