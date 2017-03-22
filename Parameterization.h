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
	HE_vert* vert;	//���ָ����Ķ���  
	HE_edge* pair;  //�ð�ߵ���һ�԰��  
	HE_face* face;  //������ڵ�ƽ��
	HE_edge* next;	//ͬ�����ϵ���һ�����
};

struct HE_vert
{
	float x, y, z;	//����ֵ
	HE_edge* edge;	//�Ӹõ�����İ��
};

struct HE_face
{
	HE_edge* edge;	//���ϵ�һ�����(��һ�����Ӧ�ı�)
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
	//	    �Բ�ͬ���͵ĵ�ִ�в���������    //
	/////////////////////////////////////

	void BorderVerticesParameterize(int TYPE, float args);		//����������߽��ϵĵ�
	void InnerVerticesParameterize(int TYPE);					//����������ڲ��ĵ�

public:

	//////////////////////
	//	 �߽�����������  //
	//////////////////////

	static const int TYPE_SQUARE = 0xFF;
	static const int TYPE_CIRCLE = 0xFE;

	//////////////////////
	//	 �ڲ������������  //
	//////////////////////

	static const int TYPE_AVERAGE_PARAMETERIZATION = 0xFD;
	static const int TYPE_SHAPE_PRESERVING_PARAMETERIZATION = 0xFC;

private:
	//void DataInitialize();									//���𽫵��������ת��Ϊ������ݽṹ
	void FindBorderInnerVertices();								//�ҵ��߽����ڲ���

private:

	////////////////////////
	//	    �������ķ���    //
	////////////////////////

	std::vector<float>* ShapePreservingParameterization(int i);

private:

	///////////////////////////////
	//    �����ȴ洢����ؾ������   //
	///////////////////////////////

	float getMatrix(float* matrix, int ld, int i, int j);
	void setMatrix(float* matrix, int ld, int i, int j, float data);

public:

	/////////////////////
	//	    ԭʼ����    //
	/////////////////////

	int m_nVertices;
	CGeoScene::CVertex* m_pVertices;
	int m_nEdges;
	CGeoScene::CEdge* m_pEdges;
	int m_nTriangles;
	CGeoScene::CTriangle* m_pTriangles;

private:

	//////////////////////
	//	    ����������    //
	//////////////////////

	std::vector<int> m_innerVertices;	
	std::vector<CPoint3D*> m_innerVerticesResult;
	std::vector<int> m_borderVertices;
	std::vector<CPoint3D*> m_borderVerticesResult;
	float m_fBorderLength;

	/////////////////////
	//	    ������    //
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