// GeneticPerceptron.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"

using namespace std;
using namespace GP;

// Функция активации нейрона. Не генерирует исключений.
double F(double _param) _NOEXCEPT
{
	return ( 1 / (1 + exp(-_param)) );
}

int _tmain(int argc, _TCHAR* argv[])
{
	const size_t COUNT_INS = 4;
	const size_t COUNT_OUTS = 3;
	const size_t COUNT_PARENTS = 100;

	// Обучим сеть на основе уроков из сгенерированного файла.
	PData lessons;
	lessons.loadFromTxtFile(COUNT_INS, COUNT_OUTS, "Weapons lessons.txt");
	lessons.printData();
	
	PWeb perceptron{ COUNT_INS, { 8, 5, COUNT_OUTS }, F };
	perceptron.bindData(&lessons);
	
	// Обучение.
	perceptron.learn(COUNT_PARENTS, 1000000, 0.01, 0, 10, 10, 100, true, 0);
	// Сохраняем сеть в файл.
	perceptron.saveToBinFile("perceptron1.pweb");
	// Загружаем сеть из файла.
	perceptron.loadFromBinFile("perceptron1.pweb");
	
	// Протестируем обученную сеть на другой выборке из той же области.
	double *helpBuffer = new double[COUNT_INS + COUNT_OUTS];
	lessons.loadFromTxtFile(COUNT_INS, COUNT_OUTS, "Weapons tests.txt");
	for (size_t e = 0; e < lessons.getCountLessons(); ++e)
	{
		for (size_t h = 0; h < COUNT_INS + COUNT_OUTS; ++h)
		{
			helpBuffer[h] = lessons.getData()[e* (COUNT_INS + COUNT_OUTS) + h];
		}
		cout << "Test " << e << ":" << endl;
		for (size_t l = 0; l < 32; ++l) cout << "-";
		cout << endl << "Correct| ";
		// Тестовые данные.
		for (size_t h = COUNT_INS; h < COUNT_INS + COUNT_OUTS; ++h)
		{
			cout.setf(ios::right);
			cout.width(5);
			cout << ((ptrdiff_t)(helpBuffer[h] * 1000)) / 1000.f << " | ";
		}
		cout << endl;
		for (size_t l = 0; l < 32; ++l) cout << "-";
		cout << endl << "  Real | ";
		perceptron.ask(helpBuffer);
		// Результаты теста.
		for (size_t h = COUNT_INS; h < COUNT_INS + COUNT_OUTS; ++h)
		{	
			cout.setf(ios::right);
			cout.width(5);
			cout << ((ptrdiff_t)(helpBuffer[h] * 1000))/1000.f << " | ";
		}
		cout << endl;
		for (size_t l = 0; l < 32; ++l) cout << "-";
		cout << endl << endl;
	}
	system("pause"); 
	return 0;
}


