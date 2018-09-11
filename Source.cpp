#if defined(UNICODE) && !defined(_UNICODE)
#define _UNICODE
#elif defined(_UNICODE) && !defined(UNICODE)
#define UNICODE
#endif
#include <tchar.h>
#include <windows.h>
#include <cmath>
#include <algorithm>
#include <stack>
#include <vector>
#include <iostream>
#include <fstream>
#include <vector>

COLORREF color = RGB(0, 0, 0);

COLORREF fillingColor = RGB(0, 0, 0);

TCHAR szClassName[] = _T("CodeBlocksWindowsApp");

struct point
{
	int x, y;
	point(int a, int b)
	{
		x = a;
		y = b;
	}
	point()
	{
		x = 0;
		y = 0;
	}
};

static point pointss[1000];
static int counterPoint = 0;
int Round(double num)
{
	return int(num + 0.5);
}

void Swap(int & f, int & s)
{
	int temp = f;
	int fNum = s;
	int sNum = temp;
}


//***********************Line************************
void DDA(HDC hdc, int xs, int ys, int xe, int ye)
{
	int dx = abs(xe - xs), dy = abs(ye - ys);
	if (dx > dy)
	{
		int xInc = (xs > xe) ? -1 : 1;
		double slope = double(dy) / double(dx), y = ys;
		double yInc = (ys > ye) ? -1 * slope : 1 * slope;
		SetPixel(hdc, xs, Round(y), color);
		while (xs != xe)
		{
			xs += xInc;
			y += yInc;
			SetPixel(hdc, xs, Round(y), color);
		}
	}
	else
	{
		int yInc = (ys > ye) ? -1 : 1;
		double slopeI = double(dx) / double(dy), x = xs;
		double xInc = (xs > xe) ? -1 * slopeI : slopeI;
		SetPixel(hdc, Round(x), ys, color);
		while (ys != ye)
		{
			x += xInc;
			ys += yInc;
			SetPixel(hdc, Round(x), ys, color);
		}
	}
}

void LineParametric(HDC hdc, int xs, int ys, int xe, int ye)
{
	int n = max(abs(xe - xs), abs(ye - ys));
	double dt = 1.0 / n;
	double dx = dt * (xe - xs), dy = dt * (ye - ys);
	double x = xs, y = ys;
	SetPixel(hdc, Round(x), Round(y), color);
	for (int i = 0; i < n; i++)
	{
		x += dx;
		y += dy;
		SetPixel(hdc, Round(x), Round(y), color);
	}
}

void LineMidPoint(HDC hdc, int xs, int xe, int ys, int ye)
{
	int dx = abs(xe - xs);
	int dy = abs(ye - ys);
	int x, y;
	bool sign = false;

	if (dy <= dx)
	{
		int d = dx - (2 * dy);
		int change_1 = -2 * dy;
		int change_2 = 2 * (dx - dy);
		if (xs > xe)
		{
			int temp = xe;
			xe = xs;
			xs = temp;

			if (ys < ye)
			{
				sign = true;
			}

			temp = ye;
			ye = ys;
			ys = temp;
		}

		if (ys > ye)
		{
			sign = true;
		}

		x = xs;
		y = ys;

		while (x < xe)
		{
			SetPixel(hdc, x, y, color);
			if (d > 0)
			{
				d += change_1;
			}
			else
			{
				d += change_2;
				if (sign == true)
				{
					y -= 1;
				}
				else
				{
					y += 1;
				}
			}
			x += 1;
		}
	}

	else
	{
		int d = dy - (2 * dx);
		int change_1 = -2 * dx;
		int change_2 = 2 * (dy - dx);
		if (ys > ye)
		{
			if (xs < xe)
			{
				sign = true;
			}
			int temp = xe;
			xe = xs;
			xs = temp;

			temp = ye;
			ye = ys;
			ys = temp;
		}

		if (xs > xe)
		{
			sign = true;
		}
		x = xs;
		y = ys;
		while (y < ye)
		{
			SetPixel(hdc, x, y, color);
			if (d > 0)
			{
				d += change_1;
			}
			else
			{
				d += change_2;
				if (sign == true)
				{
					x -= 1;
				}
				else
				{
					x += 1;
				}
			}
			y += 1;
		}
	}
}
//*****************************************************************



//***********************Circle************************
void DrawPoints(HDC hdc, int Xc, int Yc, int X, int Y, COLORREF color)
{
	SetPixel(hdc, Xc - Y, Yc + X, color);
	SetPixel(hdc, Xc + Y, Yc - X, color);
	SetPixel(hdc, Xc - Y, Yc - X, color);
	SetPixel(hdc, Xc + X, Yc + Y, color);
	SetPixel(hdc, Xc - X, Yc + Y, color);
	SetPixel(hdc, Xc + X, Yc - Y, color);
	SetPixel(hdc, Xc - X, Yc - Y, color);
	SetPixel(hdc, Xc + Y, Yc + X, color);

	pointss[counterPoint].x = Xc + Y;
	pointss[counterPoint].y = Yc - X;
	counterPoint++;
	pointss[counterPoint].x = Xc + X;
	pointss[counterPoint].y = Yc - Y;
	counterPoint++;
	pointss[counterPoint].x = Xc + X;
	pointss[counterPoint].y = Yc + Y;
	counterPoint++;
	pointss[counterPoint].x = Xc + Y;
	pointss[counterPoint].y = Yc - X;
	counterPoint++;
	pointss[counterPoint].x = Xc - Y;
	pointss[counterPoint].y = Yc + X;
	counterPoint++;
	pointss[counterPoint].x = Xc - X;
	pointss[counterPoint].y = Yc + Y;
	counterPoint++;
	pointss[counterPoint].x = Xc - Y;
	pointss[counterPoint].y = Yc - X;
	counterPoint++;
	pointss[counterPoint].x = Xc + Y;
	pointss[counterPoint].y = Yc + X;
	counterPoint++;
}
void CircleCartesian(HDC hdc, int xc, int yc, int a, int b)
{
	double rSquare = pow((double(a - xc)), 2) + pow((double(b - yc)), 2);
	int R = Round(sqrt(rSquare));
	int x = R, y = 0;
	DrawPoints(hdc, xc, yc, x, y, color);
	while (y < x)
	{
		y++;
		double p = rSquare - double(y * y);
		x = Round(sqrt(p));
		DrawPoints(hdc, xc, yc, x, y, color);
	}
}


void CirclePolar(HDC hdc, int xc, int yc, int a, int b)
{
	double rSquare = pow((double(a - xc)), 2) + pow((double(b - yc)), 2);
	int R = Round(sqrt(rSquare));
	double theta = 0, dTheta = 1.0 / R;
	int x = R, y = 0;
	DrawPoints(hdc, xc, yc, x, y, color);
	while (y < x)
	{
		theta += dTheta;
		x = Round(R * cos(theta));
		y = Round(R * sin(theta));
		DrawPoints(hdc, xc, yc, x, y, color);
	}
}


void CircleBresenham(HDC hdc, int xc, int yc, int a, int b)
{
	double rSquare = pow((double(a - xc)), 2) + pow((double(b - yc)), 2);
	int R = Round(sqrt(rSquare));
	int x = 0, y = R;
	DrawPoints(hdc, xc, yc, x, y, color);
	int decision = 1 - R;
	while (x < y)
	{
		int inside = (2 * x) + 3, outside = (2 * (x - y)) + 5;
		if (decision < 0)
		{
			decision += inside;
		}
		else
		{
			decision += outside;
			y -= 1;
		}
		x += 1;
		DrawPoints(hdc, xc, yc, x, y, color);
	}
}
//*****************************************************************






//***********************Filling************************************
int MAXENTRIES = 1000;

struct Entry
{
	int xmin, xmax;
};
void InitEntries(Entry table[])
{
	for (int i = 0; i<MAXENTRIES; i++)
	{
		table[i].xmin = MAXINT;
		table[i].xmax = -MAXINT;
	}
}
void ScanEdge(point p1, point p2, Entry table[])
{
	if (p1.y == p2.y)
		return;
	if (p1.y > p2.y)
	{
		int temp = p1.x;
		p1.x = p2.x;
		p2.x = temp;

		temp = p1.y;
		p1.y = p2.y;
		p2.y = temp;
	}
	double minv = (double)(p2.x - p1.x) / (p2.y - p1.y);
	double x = p1.x;
	int y = p1.y;
	while (y < p2.y)
	{
		if (x < table[y].xmin)
			table[y].xmin = (int)ceil(x);
		if (x  >table[y].xmax)
			table[y].xmax = (int)floor(x);
		y++;
		x += minv;
	}
}
void DrawSanLines(HDC hdc, Entry table[], COLORREF color)
{
	for (int y = 0; y < MAXENTRIES; y++)
	if (table[y].xmin < table[y].xmax)
	for (int x = table[y].xmin; x <= table[y].xmax; x++)
		SetPixel(hdc, x, y, color);
}
void ConvexFill(HDC hdc, point p[], int n, COLORREF color)
{
	Entry *table = new Entry[MAXENTRIES];
	InitEntries(table);
	point p1 = p[n - 1];
	for (int i = 0; i < n; i++)
	{
		point p2 = p[i];
		ScanEdge(p1, p2, table);
		p1 = p[i];
	}
	DrawSanLines(hdc, table, color);
	delete table;
}


template<class T>
class Stack {

private:
	T *data;
	T *curr;
	int sz;
	int items;

public:
	Stack(int s) {
		items = 0;
		data = curr = new T[sz = (s > 0 ? s : 1)];

	}


	~Stack(void) {
		if (data)
			delete[] data;

	}

	int max_stack_size(void) { return sz; }
	int stack_size(void) { return items; }
	int num(void) { return items; }

	int push(const T& a) {
		if (items < sz) {
			*curr++ = a;
			items++;

		}
		else {
			return -1;

		}
		return 0;

	}


	int dup(void) {
		if (items > 0) {
			return push(top());

		}
		else if (sz > 0) {
			curr++;
			items++;

		}
		return 0;

	}


	T& pop(void) {
		if (items > 0) {
			items--;
			return (*--curr);

		}
		else {
			return *data;

		}

	}


	T& top(void) {
		if (items > 0) {
			return *(curr - 1);

		}
		else {
			return *data;
		}
	}


};

void floodFill(HDC hdc, int x, int y, COLORREF BorderColor, COLORREF FillColor)
{
	Stack<point> st(100000);
	st.push(point(x, y));
	while (st.stack_size() != 0)
	{
		point p = st.top();
		st.pop();
		COLORREF c = GetPixel(hdc, p.x, p.y);
		if (c == BorderColor || c == FillColor)
			continue;
		SetPixel(hdc, p.x, p.y, FillColor);
		st.push(point(p.x, p.y + 1));
		st.push(point(p.x, p.y - 1));
		st.push(point(p.x - 1, p.y));
		st.push(point(p.x + 1, p.y));
	}
}
//*****************************************************************



//**************************Curves****************************
void Bezier(HDC hdc, int x1, int x2, int x3, int x4, int y1, int y2, int y3, int y4)
{
	int a1 = x1;    int a2 = 3 * (x2 - x1);    int a3 = 3 * x1 - 6 * x2 + 3 * x3;    int a4 = 3 * x2 - 3 * x3 + x4 - x1;
	int b1 = y1;    int b2 = 3 * (y2 - y1);    int b3 = 3 * y1 - 6 * y2 + 3 * y3;    int b4 = 3 * y2 - 3 * y3 + y4 - y1;
	double change = 0.001, time = 0, x, y;
	while (time < 1)
	{
		time += change;
		x = a1 + (a2*time) + (a3*pow(time, 2)) + (a4*pow(time, 3));
		y = b1 + (b2*time) + (b3*pow(time, 2)) + (b4*pow(time, 3));
		SetPixel(hdc, Round(x), Round(y), color);
	}
}

void Hermit(HDC hdc, int x1, int x2, int u1, int u2, int y1, int y2, int v1, int v2)
{
	int a1 = x1;    int a2 = u1;    int a3 = 3 * x2 - 3 * x1 - 2 * u1 - u2;    int a4 = 2 * x1 - 2 * x2 + u1 + u2;
	int b1 = y1;    int b2 = v1;    int b3 = 3 * y2 - 3 * y1 - 2 * v1 - v2;    int b4 = 2 * y1 - 2 * y2 + v1 + v2;
	double change = 0.001, time = 0, x, y;
	while (time < 1)
	{
		time += change;
		x = a1 + (a2*time) + (a3*pow(time, 2)) + (a4*pow(time, 3));
		y = b1 + (b2*time) + (b3*pow(time, 2)) + (b4*pow(time, 3));
		SetPixel(hdc, Round(x), Round(y), color);
	}
}
//*****************************************************************




//***********************Clipping****************************
union OutCode
{
	unsigned All : 4;
	struct{ unsigned left : 1, top : 1, right : 1, bottom : 1; };
};
OutCode GetOutCode(double x, double y, int xleft, int ytop, int xright, int ybottom)
{
	OutCode out;
	out.All = 0;
	if (x<xleft)out.left = 1; else if (x>xright)out.right = 1;
	if (y<ytop)out.top = 1; else if (y>ybottom)out.bottom = 1;
	return out;
}
void VIntersect(double xs, double ys, double xe, double ye, int x, double *xi, double *yi)
{
	*xi = x;
	*yi = ys + (x - xs)*(ye - ys) / (xe - xs);
}
void HIntersect(double xs, double ys, double xe, double ye, int y, double *xi, double *yi)
{
	*yi = y;
	*xi = xs + (y - ys)*(xe - xs) / (ye - ys);
}
void CohenSuth(HDC hdc, int xs, int ys, int xe, int ye, int xleft, int ytop, int xright, int ybottom)
{
	double x1 = xs, y1 = ys, x2 = xe, y2 = ye;
	OutCode out1 = GetOutCode(x1, y1, xleft, ytop, xright, ybottom);
	OutCode out2 = GetOutCode(x2, y2, xleft, ytop, xright, ybottom);
	while ((out1.All || out2.All) && !(out1.All & out2.All))
	{
		double xi, yi;
		if (out1.All)
		{
			if (out1.left)VIntersect(x1, y1, x2, y2, xleft, &xi, &yi);
			else if (out1.top)HIntersect(x1, y1, x2, y2, ytop, &xi, &yi);
			else if (out1.right)VIntersect(x1, y1, x2, y2, xright, &xi, &yi);
			else HIntersect(x1, y1, x2, y2, ybottom, &xi, &yi);
			x1 = xi;
			y1 = yi;
			out1 = GetOutCode(x1, y1, xleft, ytop, xright, ybottom);
		}
		else
		{
			if (out2.left)VIntersect(x1, y1, x2, y2, xleft, &xi, &yi);
			else if (out2.top)HIntersect(x1, y1, x2, y2, ytop, &xi, &yi);
			else if (out2.right)VIntersect(x1, y1, x2, y2, xright, &xi, &yi);
			else HIntersect(x1, y1, x2, y2, ybottom, &xi, &yi);
			x2 = xi;
			y2 = yi;
			out2 = GetOutCode(x2, y2, xleft, ytop, xright, ybottom);
		}
	}
	if (!out1.All && !out2.All)
	{
		MoveToEx(hdc, Round(x1), Round(y1), NULL);
		LineTo(hdc, Round(x2), Round(y2));
	}
}



struct Vertex
{
	double x, y;
	Vertex(int x1 = 0, int y1 = 0)
	{
		x = x1;
		y = y1;
	}
};

typedef std::vector<Vertex> VertexList;
typedef bool(*IsInFunc)(Vertex& v, int edge);
typedef Vertex(*IntersectFunc)(Vertex& v1, Vertex& v2, int edge);

bool InLeft(Vertex& v, int edge)
{
	return v.x >= edge;
}
bool InRight(Vertex& v, int edge)
{
	return v.x <= edge;
}
bool InTop(Vertex& v, int edge)
{
	return v.y >= edge;
}
bool InBottom(Vertex& v, int edge)
{
	return v.y <= edge;
}
Vertex VIntersect(Vertex& v1, Vertex& v2, int xedge)
{
	Vertex res;
	res.x = xedge;
	res.y = v1.y + (xedge - v1.x)*(v2.y - v1.y) / (v2.x - v1.x);
	return res;
}
Vertex HIntersect(Vertex& v1, Vertex& v2, int yedge)
{
	Vertex res;
	res.y = yedge;
	res.x = v1.x + (yedge - v1.y)*(v2.x - v1.x) / (v2.y - v1.y);
	return res;
}
VertexList ClipWithEdge(VertexList p, int edge, IsInFunc In, IntersectFunc Intersect)
{
	VertexList OutList;
	Vertex v1 = p[p.size() - 1];
	bool v1_in = In(v1, edge);
	for (int i = 0; i<(int)p.size(); i++)
	{
		Vertex v2 = p[i];
		bool v2_in = In(v2, edge);
		if (!v1_in && v2_in)
		{
			OutList.push_back(Intersect(v1, v2, edge));
			OutList.push_back(v2);
		}
		else if (v1_in && v2_in) OutList.push_back(v2);
		else if (v1_in) OutList.push_back(Intersect(v1, v2, edge));
		v1 = v2;
		v1_in = v2_in;
	}
	return OutList;
}

void PolygonClip(HDC hdc, point p[], int n, int xleft, int ytop, int xright, int ybottom)
{
	VertexList vlist;
	for (int i = 0; i<n; i++)vlist.push_back(Vertex(p[i].x, p[i].y));
	vlist = ClipWithEdge(vlist, xleft, InLeft, VIntersect);
	vlist = ClipWithEdge(vlist, ytop, InTop, HIntersect);
	vlist = ClipWithEdge(vlist, xright, InRight, VIntersect);
	vlist = ClipWithEdge(vlist, ybottom, InBottom, HIntersect);
	Vertex v1 = vlist[vlist.size() - 1];
	for (int i = 0; i<(int)vlist.size(); i++)
	{
		Vertex v2 = vlist[i];
		MoveToEx(hdc, Round(v1.x), Round(v1.y), NULL);
		LineTo(hdc, Round(v2.x), Round(v2.y));
		v1 = v2;
	}
}


//*****************************************************************

//**************************Save And Load******************************
typedef unsigned short     uint16_t;
typedef unsigned int       uint32_t;

bool HDCToFile(const char* FilePath, HDC Context, RECT Area, uint16_t BitsPerPixel = 24)
{
	uint32_t Width = Area.right - Area.left;
	uint32_t Height = Area.bottom - Area.top;
	BITMAPINFO Info;
	BITMAPFILEHEADER Header;
	memset(&Info, 0, sizeof(Info));
	memset(&Header, 0, sizeof(Header));
	Info.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
	Info.bmiHeader.biWidth = Width;
	Info.bmiHeader.biHeight = Height;
	Info.bmiHeader.biPlanes = 1;
	Info.bmiHeader.biBitCount = BitsPerPixel;
	Info.bmiHeader.biCompression = BI_RGB;
	Info.bmiHeader.biSizeImage = Width * Height * (BitsPerPixel > 24 ? 4 : 3);
	Header.bfType = 0x4D42;
	Header.bfOffBits = sizeof(BITMAPFILEHEADER)+sizeof(BITMAPINFOHEADER);
	char* Pixels = NULL;
	HDC MemDC = CreateCompatibleDC(Context);
	HBITMAP Section = CreateDIBSection(Context, &Info, DIB_RGB_COLORS, (void**)&Pixels, 0, 0);
	BITMAP qBitmap;
	int iReturn = GetObject(reinterpret_cast<HGDIOBJ>(Section), sizeof(BITMAP),
		reinterpret_cast<LPVOID>(&qBitmap));
	DeleteObject(SelectObject(MemDC, Section));
	BitBlt(MemDC, 0, 0, qBitmap.bmWidth, qBitmap.bmHeight, Context, 0, 0, SRCCOPY);
	DeleteDC(MemDC);
	std::fstream hFile(FilePath, std::ios::out | std::ios::binary);
	if (hFile.is_open())
	{
		hFile.write((char*)&Header, sizeof(Header));
		hFile.write((char*)&Info.bmiHeader, sizeof(Info.bmiHeader));
		hFile.write(Pixels, (((BitsPerPixel * Width + 31) & ~31) / 8) * Height);
		hFile.close();
		DeleteObject(Section);
		return true;
	}
	DeleteObject(Section);
	return false;
}

bool LoadAndBlitBitmap(LPCSTR szFileName, HDC hWinDC)
{
	// Load the bitmap image file
	HBITMAP hBitmap;
	hBitmap = (HBITMAP)::LoadImage(NULL, (LPCWSTR)szFileName, IMAGE_BITMAP, 0, 0,
		LR_LOADFROMFILE);
	// Verify that the image was loaded
	if (hBitmap == NULL) {
		::MessageBox(NULL, (LPCWSTR)L"LoadImage Failed", (LPCWSTR)L"Error", MB_OK);
		return false;
	}

	// Create a device context that is compatible with the window
	HDC hLocalDC;
	hLocalDC = ::CreateCompatibleDC(hWinDC);
	// Verify that the device context was created
	if (hLocalDC == NULL) {
		::MessageBox(NULL, (LPCWSTR)("CreateCompatibleDC Failed"), (LPCWSTR)("Error"), MB_OK);
		return false;
	}

	// Get the bitmap's parameters and verify the get
	BITMAP qBitmap;
	int iReturn = GetObject(reinterpret_cast<HGDIOBJ>(hBitmap), sizeof(BITMAP),
		reinterpret_cast<LPVOID>(&qBitmap));
	if (!iReturn) {
		::MessageBox(NULL, (LPCWSTR)L"GetObject Failed", (LPCWSTR)L"Error", MB_OK);
		return false;
	}

	// Select the loaded bitmap into the device context
	HBITMAP hOldBmp = (HBITMAP)::SelectObject(hLocalDC, hBitmap);
	if (hOldBmp == NULL) {
		::MessageBox(NULL, (LPCWSTR)L"SelectObject Failed", (LPCWSTR)L"Error", MB_OK);
		return false;
	}

	// Blit the dc which holds the bitmap onto the window's dc
	BOOL qRetBlit = ::BitBlt(hWinDC, 0, 0, qBitmap.bmWidth, qBitmap.bmHeight,
		hLocalDC, 0, 0, SRCCOPY);
	if (!qRetBlit) {
		::MessageBox(NULL, (LPCWSTR)L"Blit Failed", (LPCWSTR)L"Error", MB_OK);
		return false;
	}

	// Unitialize and deallocate resources
	::SelectObject(hLocalDC, hOldBmp);
	::DeleteDC(hLocalDC);
	::DeleteObject(hBitmap);
	return true;
}

void save(HWND hWnd)
{
	HDC hdc = GetDC(hWnd);
	RECT rect;
	GetWindowRect(hWnd, &rect);
	HDCToFile("save.bmp", hdc, rect);
}
//*****************************************************************


HMENU menue = CreateMenu();
HMENU lineMenue = CreateMenu();
HMENU circleMenue = CreateMenu();
HMENU fillingMenue = CreateMenu();
HMENU curveMenue = CreateMenu();
HMENU clippingMenue = CreateMenu();
HMENU fileMenue = CreateMenu();
HMENU coloringMenue = CreateMenu();
HMENU coloringFillingMenue = CreateMenu();
int xs, ys, xe, ye, R, checkShape, counter = 0, bezierOrHermitCounter = 0;
int xx1, xx2, xx3, xx4, yy1, yy2, yy3, yy4;
int xleft, ytop, xright, ybottom;
LRESULT WINAPI MyWindowProc(HWND hWnd, UINT m, WPARAM wp, LPARAM lp)
{
	HDC hdc;
	hdc = GetDC(hWnd);
	
	switch (m)
	{
	case WM_CREATE:
		//file
		AppendMenu(menue, MF_POPUP, (UINT_PTR)fileMenue, L"File");
		AppendMenu(fileMenue, MF_STRING, 1, L"Save");
		AppendMenu(fileMenue, MF_STRING, 2, L"Load");

		//line
		AppendMenu(menue, MF_POPUP, (UINT_PTR)lineMenue, L"Line");
		AppendMenu(lineMenue, MF_STRING, 3, L"DDA");
		AppendMenu(lineMenue, MF_STRING, 4, L"Mid Point");
		AppendMenu(lineMenue, MF_STRING, 5, L"Parametric");


		//circle
		AppendMenu(menue, MF_POPUP, (UINT_PTR)circleMenue, L"Circle");
		AppendMenu(circleMenue, MF_STRING, 6, L"Cartesian");
		AppendMenu(circleMenue, MF_STRING, 7, L"Polar");
		AppendMenu(circleMenue, MF_STRING, 8, L"Mid Point");


		//filling
		AppendMenu(menue, MF_POPUP, (UINT_PTR)fillingMenue, L"Filling");
		AppendMenu(fillingMenue, MF_STRING, 9, L"Convex");
		AppendMenu(fillingMenue, MF_STRING, 10, L"Fill Shape");


		//curve
		AppendMenu(menue, MF_POPUP, (UINT_PTR)curveMenue, L"Curve");
		AppendMenu(curveMenue, MF_STRING, 11, L"Bezeir");
		AppendMenu(curveMenue, MF_STRING, 12, L"Hermit");


		//Clipping
		AppendMenu(menue, MF_POPUP, (UINT_PTR)clippingMenue, L"Clipping");
		AppendMenu(clippingMenue, MF_STRING, 13, L"Line");
		AppendMenu(lineMenue, MF_STRING, 14, NULL);
		AppendMenu(clippingMenue, MF_STRING, 15, L"Polygon");

		
		AppendMenu(fillingMenue, MF_STRING, 16, L"Convex Circle");
		AppendMenu(clippingMenue, MF_STRING, 17, L"Shape");

		AppendMenu(menue, MF_POPUP, (UINT_PTR)coloringMenue, L"Colors");
		AppendMenu(coloringMenue, MF_STRING, 18, L"White");
		AppendMenu(coloringMenue, MF_STRING, 19, L"Red");
		AppendMenu(coloringMenue, MF_STRING, 20, L"Black");

		AppendMenu(menue, MF_POPUP, (UINT_PTR)coloringFillingMenue, L"Filling Colors");
		AppendMenu(coloringFillingMenue, MF_STRING, 21, L"White");
		AppendMenu(coloringFillingMenue, MF_STRING, 22, L"Red");
		AppendMenu(coloringFillingMenue, MF_STRING, 23, L"Black");
		AppendMenu(coloringFillingMenue, MF_STRING, 24, L"Yellow");

		SetMenu(hWnd, menue);
		break;
	case WM_COMMAND:
		if (wp == 1)
		{
			save(hWnd);
		}
		if (wp == 2)
		{
			checkShape = 2;
			LoadAndBlitBitmap((LPCSTR)L"save.bmp", hdc);
		}
			
		if (wp == 3)
			checkShape = 3;
		if (wp == 4)
			checkShape = 4;
		if (wp == 5)
			checkShape = 5;
		if (wp == 6)
			checkShape = 6;
		if (wp == 7)
			checkShape = 7;
		if (wp == 8)
			checkShape = 8;
		if (wp == 9)
		{
			checkShape = 9;
			counterPoint = 0;
			memset(pointss, 0, sizeof(point));
		}
		if (wp == 10)
			checkShape = 10;
		if (wp == 11)
			checkShape = 11;
		if (wp == 12)
			checkShape = 12;
		if (wp == 13)
			checkShape = 13;
		if (wp == 14)
			checkShape = 14;
		if (wp == 15)
		{
			checkShape = 15;
			counterPoint = 0;
			memset(pointss, 0, sizeof(point));
		}
		if (wp == 16)
		{
			ConvexFill(hdc, pointss, counterPoint, fillingColor);
			counterPoint = 0;
			memset(pointss, 0, sizeof(point));
		}
		if (wp == 17)
		{
			checkShape = 17;
		}
		if (wp == 18)
		{
			color = RGB(255, 255, 255);
		}
		if (wp == 19)
		{
			color = RGB(255, 0, 0);
		}
		if (wp == 20)
		{
			color = RGB(0, 0, 0);
		}
		if (wp == 21)
		{
			fillingColor = RGB(255, 255, 255);
		}
		if (wp == 22)
		{
			fillingColor = RGB(255, 0, 0);
		}
		if (wp == 23)
		{
			fillingColor = RGB(232, 245, 17);
		}
		if (wp == 24)
		{
			fillingColor = RGB(245, 39, 17);
		}
		break;
	case WM_LBUTTONDOWN:
		counter++;
		if (counter % 2 != 0)
		{
			xs = LOWORD(lp);
			ys = HIWORD(lp);
			pointss[counterPoint].x = LOWORD(lp);
			pointss[counterPoint].y = HIWORD(lp);
			counterPoint++;
			if (checkShape == 11)
			{
				bezierOrHermitCounter++;
				if (bezierOrHermitCounter == 1)
				{
					xx1 = LOWORD(lp);
					yy1 = HIWORD(lp);
				}
				else if (bezierOrHermitCounter == 3)
				{
					xx3 = LOWORD(lp);
					yy3 = HIWORD(lp);
				}
			}
			if (checkShape == 12)
			{
				bezierOrHermitCounter++;
				if (bezierOrHermitCounter == 1)
				{
					xx1 = LOWORD(lp);
					yy1 = HIWORD(lp);
				}
				else if (bezierOrHermitCounter == 3)
				{
					xx3 = LOWORD(lp);
					yy3 = HIWORD(lp);
				}
			}
			if (checkShape == 13)
			{
				xleft = LOWORD(lp);
				ytop = HIWORD(lp);
			}
			if (checkShape == 15)
			{
				xleft = LOWORD(lp);
				ytop = HIWORD(lp);
				counterPoint = 0;
			}
			else if (checkShape == 18)
			{
				xleft = LOWORD(lp);
				ytop = HIWORD(lp);
			}
		}
		else if (counter % 2 == 0)
		{
			xe = LOWORD(lp);
			ye = HIWORD(lp);

			if (checkShape == 3)
			{
				DDA(hdc, xs, ys, xe, ye);
			}
			else if (checkShape == 4)
			{
				LineMidPoint(hdc, xs, xe, ys, ye);
			}
			else if (checkShape == 5)
			{
				LineParametric(hdc, xs, ys, xe, ye);
			}
			else if (checkShape == 6)
			{
				CircleCartesian(hdc, xs, ys, xe, ye);
			}
			else if (checkShape == 7)
			{
				CirclePolar(hdc, xs, ys, xe, ye);
			}
			else if (checkShape == 8)
			{
				CircleBresenham(hdc, xs, ys, xe, ye);
			}
			else if (checkShape == 9)
			{
				pointss[counterPoint].x = LOWORD(lp);
				pointss[counterPoint].y = HIWORD(lp);
				counterPoint++;
			}
			else if (checkShape == 10)
			{
				floodFill(hdc, xs, ys, color, fillingColor);
			}
			else if (checkShape == 11)
			{
				bezierOrHermitCounter++;
				if (bezierOrHermitCounter == 2)
				{
					xx2 = LOWORD(lp);
					yy2 = HIWORD(lp);
				}
				if (bezierOrHermitCounter == 4)
				{
					xx4 = LOWORD(lp);
					yy4 = HIWORD(lp);
					Bezier(hdc, xx1, xx2, xx3, xx4, yy1, yy2, yy3, yy4);
					bezierOrHermitCounter = 0;
				}

			}
			else if (checkShape == 12)
			{
				bezierOrHermitCounter++;
				if (bezierOrHermitCounter == 2)
				{
					xx2 = LOWORD(lp);
					yy2 = HIWORD(lp);
				}
				if (bezierOrHermitCounter == 4)
				{
					xx4 = LOWORD(lp);
					yy4 = HIWORD(lp);
					Hermit(hdc, xx1, xx2, xx3, xx4, yy1, yy2, yy3, yy4);
					bezierOrHermitCounter = 0;
				}

			}
			else if (checkShape == 13)
			{
				xright = LOWORD(lp);
				ybottom = HIWORD(lp);
				Rectangle(hdc, xleft, ytop, xright, ybottom);
				checkShape = 14;
			}
			else if (checkShape == 14)
			{
				CohenSuth(hdc, xs, ys, xe, ye, xleft, ytop, xright, ybottom);
			}
			else if (checkShape == 15)
			{
				xright = LOWORD(lp);
				ybottom = HIWORD(lp);
				Rectangle(hdc, xleft, ytop, xright, ybottom);;
				counterPoint = 0;
				checkShape = 16;
			}
			else if (checkShape == 16)
			{
				pointss[counterPoint].x = LOWORD(lp);
				pointss[counterPoint].y = HIWORD(lp);
				counterPoint++;
			}
			else if (checkShape == 17)
			{
				pointss[counterPoint].x = LOWORD(lp);
				pointss[counterPoint].y = HIWORD(lp);
				counterPoint++;
			}
			else if (checkShape == 18)
			{
				xright = LOWORD(lp);
				ybottom = HIWORD(lp);
				Rectangle(hdc, xleft, ytop, xright, ybottom);;
				PolygonClip(hdc, pointss, counterPoint, xleft, ytop, xright, ybottom);
				counterPoint = 0;
				memset(pointss, 0, sizeof(point));
			}
			InvalidateRect(hWnd, NULL, FALSE);
		}
		break;
	case WM_RBUTTONDOWN:
		if (checkShape == 9)
		{
			ConvexFill(hdc, pointss, counterPoint, RGB(0, 0, 0));
			counterPoint = 0;
			memset(pointss, 0, sizeof(point));
		}
		if (checkShape == 16)
		{
			PolygonClip(hdc, pointss, counterPoint, xleft, ytop, xright, ybottom);
			counterPoint = 0;
			memset(pointss, 0, sizeof(point));
		}
		if (checkShape == 17)
		{
			checkShape = 18;
		}
		
		break;
	case WM_DESTROY:
		PostQuitMessage(0);
		break;
	default: return DefWindowProc(hWnd, m, wp, lp);
	}
	return 0;

}
int WINAPI WinMain(HINSTANCE hThisInstance, HINSTANCE hPrevInstance, LPSTR lpszArgument, int nCmdShow)
{
	HWND hwnd;
	MSG messages;
	WNDCLASSEX wincl;
	wincl.hInstance = hThisInstance;
	wincl.lpszClassName = szClassName;
	wincl.lpfnWndProc = MyWindowProc;
	wincl.style = CS_DBLCLKS;
	wincl.cbSize = sizeof (WNDCLASSEX);
	wincl.hIcon = LoadIcon(NULL, IDI_APPLICATION);
	wincl.hIconSm = LoadIcon(NULL, IDI_APPLICATION);
	wincl.hCursor = LoadCursor(NULL, IDC_ARROW);
	wincl.lpszMenuName = NULL;
	wincl.cbClsExtra = 0;
	wincl.cbWndExtra = 0;
	wincl.hbrBackground = (HBRUSH)COLOR_BACKGROUND;
	if (!RegisterClassEx(&wincl))
		return 0;
	hwnd = CreateWindowEx(0, szClassName, _T("Graphics Algorithms"), WS_OVERLAPPEDWINDOW, CW_USEDEFAULT, CW_USEDEFAULT, 644, 475, HWND_DESKTOP, NULL, hThisInstance, NULL);
	ShowWindow(hwnd, nCmdShow);
	while (GetMessage(&messages, NULL, 0, 0))
	{
		TranslateMessage(&messages);
		DispatchMessage(&messages);
	}
	return messages.wParam;
}