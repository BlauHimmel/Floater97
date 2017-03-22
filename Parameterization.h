#ifndef PARAMETERIZATION
#define PARAMETERIZATION

#include "GeoScene.h"
#include <vector>
#include <gl\glut.h>
#include <gl\GL.h>
#include <gl\GLU.h>
#include <gl\GLAUX.H>

/*
struct HE_edge
{	
	HE_vert* vert;	//半边指向方向的顶点  
	HE_edge* pair;  //该半边的另一对半边  
	HE_face* face;  //半边所在的平面
	HE_edge* next;	//同方向上的下一个半边
};

struct HE_vert
{
	float x, y, z;	//坐标值
	HE_edge* edge;	//从该点射出的半边
};

struct HE_face
{
	HE_edge* edge;	//面上的一条半边(第一个点对应的边)
};
*/

#define PI 3.1415926

class CParameterization
{
public:
	CParameterization(){};
	CParameterization(int nVertices, int nEdges, int nTriangles,
		CGeoScene::CVertex* pVertices, CGeoScene::CEdge* pEdges,
		CGeoScene::CTriangle* pTriangles) :m_nVertices(nVertices), 
		m_nEdges(nEdges), m_nTriangles(nTriangles)
		{
			m_pVertices = pVertices;
			m_pEdges = pEdges;
			m_pTriangles = pTriangles;
			m_pVerticesTrans = nullptr;
		}
	CParameterization(CParameterization& other);
	~CParameterization();

public:
	CPoint3D** Calculate(float args);
    void Draw(float args);
	void Draw(float args, const char* filename);
	void loadImage(const char *filename);

private:
	GLbyte* gltloadbmp(const char *filename, GLint *width, GLint *height);

public:
	static void* __pImage;
	static GLint __width;
	static GLint __height;

public:

	/////////////////////////////////////
	//	    对不同类型的点执行参数化方法    //
	/////////////////////////////////////

	void BorderVerticesParameterize(int TYPE, float args);		//负责参数化边界上的点
	void InnerVerticesParameterize(int TYPE);					//负责参数化内部的点

public:

	//////////////////////
	//	 边界点参数化类型  //
	//////////////////////

	static const int TYPE_SQUARE = 0xFF;
	static const int TYPE_CIRCLE = 0xFE;

	//////////////////////
	//	 内部点参数化类型  //
	//////////////////////

	static const int TYPE_AVERAGE_PARAMETERIZATION = 0xFD;
	static const int TYPE_SHAPE_PRESERVING_PARAMETERIZATION = 0xFC;

private:
	//void DataInitialize();									//负责将导入的数据转换为半边数据结构
	void FindBorderInnerVertices();								//找到边界点和内部点

private:

	////////////////////////
	//	    参数化的方法    //
	////////////////////////

	std::vector<float>* ShapePreservingParameterization(int i);

private:

	///////////////////////////////
	//    行优先存储的相关矩阵操作   //
	///////////////////////////////

	float getMatrix(float* matrix, int ld, int i, int j);
	void setMatrix(float* matrix, int ld, int i, int j, float data);

public:

	/////////////////////
	//	    原始数据    //
	/////////////////////

	int m_nVertices;
	CGeoScene::CVertex* m_pVertices;
	int m_nEdges;
	CGeoScene::CEdge* m_pEdges;
	int m_nTriangles;
	CGeoScene::CTriangle* m_pTriangles;

private:

	//////////////////////
	//	    参数化数据    //
	//////////////////////

	std::vector<int> m_innerVertices;	
	std::vector<CPoint3D*> m_innerVerticesResult;
	std::vector<int> m_borderVertices;
	std::vector<CPoint3D*> m_borderVerticesResult;
	float m_fBorderLength;

	/////////////////////
	//	    处理结果    //
	/////////////////////
	CPoint3D** m_pVerticesTrans;

/*
private:
	std::vector<HE_vert*> innerVertices;
	std::vector<HE_vert*> borderVertices;
	std::vector<HE_vert*> vertices;
	std::vector<HE_edge*> edges;
	std::vector<HE_face*> faces;
*/
};

#endif