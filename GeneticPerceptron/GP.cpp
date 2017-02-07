#include "stdafx.h"

using namespace std;
using namespace GP;

// Конструктор с параметрами. Генерирует исключение, если: 
// -_countWeights < 1.
PNeuron::PNeuron(const size_t _countWeights) : countWeights{ _countWeights }, weights{ nullptr }
{
	if (countWeights < 1)
	{
		throw(std::exception("_countWeights < 1"));
	}
	else {
		weights = new double[countWeights];
		for (size_t w = 0; w < countWeights; ++w)
		{
			weights[w] = 0;
		}
	}
}
// Конструктор копирования. Не генерирует исключений.
PNeuron::PNeuron(const PNeuron &_object) : countWeights{ _object.countWeights }, weights{ nullptr }
{
	if (countWeights != 0)
	{
		weights = new double[countWeights];
		memcpy(weights, _object.weights, sizeof(double)* countWeights);
	}
}
// Конструктор перемещения. Не генерирует исключений.
PNeuron::PNeuron(PNeuron &&_object) _NOEXCEPT : countWeights{ _object.countWeights }, weights{ _object.weights }
{
	_object.countWeights = 0;
	_object.weights = nullptr;
}
// Копирующий оператор присваивания, полностью оптимизированный, перестраивает топологию только при необходимости. Не генерирует исключений.
PNeuron &PNeuron::operator=(const PNeuron &_object) _NOEXCEPT
{
	if (this != &_object)
	{
		if (countWeights == _object.countWeights)
		{
			if (countWeights == 0)
			{
				weights = nullptr;
			}
			else {
				memcpy(weights, _object.weights, sizeof(double)* countWeights);
			}
		}
		else {
			delete[] weights;
			countWeights = _object.countWeights;
			if (countWeights == 0)
			{
				weights = nullptr;
			}
			else {
				weights = new double[countWeights];
				memcpy(weights, _object.weights, sizeof(double)* countWeights);
			}
		}
	}
	return *this;
}
// Перемещающий оператор присваивания. Не генерирует исключений.
PNeuron &PNeuron::operator=(PNeuron &&_object) _NOEXCEPT
{
	if (this != &_object)
	{
		if (countWeights == 0)
		{
			delete[] weights;
		}
		countWeights = _object.countWeights;
		weights = _object.weights;
		_object.countWeights = 0;
		_object.weights = nullptr;
	}
	return *this;
}
// Деструктор с полным занулением . Не генерирует исключений.
PNeuron::~PNeuron() _NOEXCEPT
{
	countWeights = 0;
	delete[] weights;
	weights = nullptr;
}
// Пропускает данные через нейрон, помещая сгенерированное значение в _out, без проверки размерности входных данных. Не генерирует исключений.
void PNeuron::genOut(const double *_ins, double &_out, double(*_f)(double)) _NOEXCEPT
{
	_out = 0;
	for (size_t w = 0; w < countWeights; ++w)
	{
		_out += _ins[w] * weights[w];
	}
	if (_f == nullptr)
	{
		_out = _out / countWeights;
	}
	else {
		_out = _f(_out / countWeights);
	}
}
// Возвращает ссылку на вес с индексом _w, не проверяя границы массива. Не генерирует исключений.
double &PNeuron::operator[](size_t _w) _NOEXCEPT
{
	return weights[_w];
}
// Возвращает ссылку на число связей в нейроне. Число связей нельзя изменить. Не генерирует исключений.
const size_t &PNeuron::getCountWeights() const _NOEXCEPT
{
	return countWeights;
}
// Конструктор с параметрами. Во время создания, все веса всех нейронов сети заполняются 0. 
// Генерирует исключение, если: 
// - _topology.size() < 1;
// - _topology[t] < 1;
// - _countIns < 1.
PWeb::PWeb(const size_t _countIns, initializer_list<size_t> _topology, double(*_f)(double)) : countIns{ _countIns },
countLayers{ 0 }, topology{ nullptr }, neurons{ nullptr }, sumSigma{ 0.0 }, f{ _f }, data{ nullptr }
{
	if (_topology.size() < 1) throw(std::exception("_topology.size() < 1"));
	for (auto &t : _topology)
	{
		if (t < 1) throw(std::exception("_topology[t] < 1"));
	}
	if (_countIns < 1) throw(std::exception("_countIns < 1"));
	countLayers = _topology.size();
	topology = new size_t[countLayers];
	// helpIns и helpOuts имеют одинаковый размер = countHelpValues. Это нужно для использования быстрого std::swap(helpIns, helpOuts), 
	// вместо memcpy(). Небольшой перерасход памяти, но зато значительный выигрыш в скорости в методах PWeb::test() и PWeb::ask().
	countHelpValues = _countIns;
	size_t i = 0;
	for (auto &t : _topology)
	{
		topology[i] = t;
		++i;
		if (countHelpValues < t) countHelpValues = t;
	}
	helpIns = new double[countHelpValues];
	helpOuts = new double[countHelpValues];
	neurons = new PNeuron**[countLayers];
	for (size_t l = 0; l < countLayers; ++l)
	{
		neurons[l] = new PNeuron*[topology[l]];
		for (size_t n = 0; n < topology[l]; ++n)
		{
			if (l == 0)
			{
				neurons[l][n] = new PNeuron(countIns);
			}
			else {
				neurons[l][n] = new PNeuron(topology[l - 1]);
			}
		}
	}
}
// Конструктор копирования. Не генерирует исключений. 
PWeb::PWeb(const PWeb &_object) _NOEXCEPT : f{ _object.f }, countIns{ _object.countIns }, 
countLayers{ _object.countLayers }, topology{ nullptr }, countHelpValues{ _object.countHelpValues }, neurons{ nullptr }, 
sumSigma{ _object.sumSigma }, data{ _object.data }
{
	topology = new size_t[countLayers];
	memcpy(topology, _object.topology, sizeof(size_t) * countLayers);
	helpIns = new double[countHelpValues];
	helpOuts = new double[countHelpValues];
	neurons = new PNeuron**[countLayers];
	for (size_t l = 0; l < countLayers; ++l)
	{
		neurons[l] = new PNeuron*[topology[l]];
		for (size_t n = 0; n < topology[l]; ++n)
		{
			neurons[l][n] = new PNeuron(*_object.neurons[l][n]);
		}
	}
}
// Конструктор перемещения. Не генерирует исключений.
PWeb::PWeb(PWeb &&_object) _NOEXCEPT : countIns{ _object.countIns }, countLayers{ _object.countLayers }, topology{ _object.topology },
countHelpValues{ _object.countHelpValues }, helpIns{ _object.helpIns }, helpOuts{ _object.helpOuts }, neurons{ _object.neurons }, 
sumSigma{ _object.sumSigma }, f{ _object.f }, data{ _object.data }
{
	_object.countIns = 0;
	_object.countLayers = 0;
	_object.topology = nullptr;
	_object.countHelpValues = 0;
	_object.helpIns = nullptr;
	_object.helpOuts = nullptr;
	_object.neurons = nullptr;
	_object.sumSigma = 0.0;
	_object.f = nullptr;
	_object.data = nullptr;
}
// Копирующий оператор присваивания, полностью оптимизирован, перестраивает топологию только при необходимости, 
// но перестраивает полностью. Не генерирует исключений.
PWeb &PWeb::operator=(const PWeb &_object) _NOEXCEPT
{
	if (this != &_object)
	{
		bool equal = true;
		if (countLayers == _object.countLayers)
		{
			for (size_t l = 0; l < countLayers; ++l)
			{
				if (topology[l] != _object.topology[l])
				{
					equal = false;
					goto A;
				}
			}
		}
		else {
			equal = false;
		}
	A:
		if (equal == true)
		{
			// Топология совпала:
			for (size_t l = 0; l < countLayers; ++l)
			{
				for (size_t n = 0; n < topology[l]; ++n)
				{
					*neurons[l][n] = *_object.neurons[l][n];
				}
			}
		}
		else {
			// Топология не совпала:
			for (size_t l = 0; l < countLayers; ++l)
			{
				for (size_t n = 0; n < topology[l]; ++n)
				{
					delete neurons[l][n];
				}
				delete[] neurons[l];
			}
			delete neurons;
			delete[] topology;
			delete[] helpIns;
			delete[] helpOuts;
			countIns = _object.countIns;
			countLayers = _object.countLayers;
			topology = new size_t[countLayers];
			memcpy(topology, _object.topology, sizeof(size_t) * countLayers);
			countHelpValues = _object.countHelpValues;
			helpIns = new double[countHelpValues];
			helpOuts = new double[countHelpValues];
			neurons = new PNeuron**[countLayers];
			for (size_t l = 0; l < countLayers; ++l)
			{
				neurons[l] = new PNeuron*[topology[l]];
				for (size_t n = 0; n < topology[l]; ++n)
				{
					neurons[l][n] = new PNeuron(*_object.neurons[l][n]);
				}
			}
			
		}
		sumSigma = _object.sumSigma;
		f = _object.f;
		data = _object.data;
	}
	return *this;
}
// Перемещающий оператор присваивания. Не генерирует исключений.
PWeb &PWeb::operator=(PWeb &&_object) _NOEXCEPT
{
	if (this != &_object)
	{
		for (size_t l = 0; l < countLayers; ++l)
		{
			for (size_t n = 0; n < topology[l]; ++n)
			{
				delete neurons[l][n];
			}
			delete[] neurons[l];
		}
		delete[] neurons;
		delete[] topology;
		delete[] helpIns;
		delete[] helpOuts;
		countIns = _object.countIns;
		countLayers = _object.countLayers;
		topology = _object.topology;
		countHelpValues = _object.countHelpValues;
		helpIns = _object.helpIns;
		helpOuts = _object.helpOuts;
		neurons = _object.neurons;
		sumSigma = _object.sumSigma;
		f = _object.f;
		data = _object.data;
		_object.countIns = 0;
		_object.countLayers = 0;
		_object.topology = nullptr;
		_object.countHelpValues = 0;
		_object.helpIns = nullptr;
		_object.helpOuts = nullptr;
		_object.neurons = nullptr;
		_object.sumSigma = 0.0;	
		_object.f = nullptr;
		_object.data = nullptr;
	}
	return *this;
}
// Деструктор с полным занулением. Не генерирует исключений.
PWeb::~PWeb() _NOEXCEPT
{
	f = nullptr;
	sumSigma = +0.0;
	for (size_t l = 0; l < countLayers; ++l)
	{
		for (size_t n = 0; n < topology[l]; ++n)
		{
			delete neurons[l][n];
			neurons[l][n] = nullptr;
		}
		delete[] neurons[l];
		neurons[l] = nullptr;
	}
	delete[] neurons;
	neurons = nullptr;
	delete[] topology;
	topology = nullptr;
	delete[] helpOuts;
	helpOuts = nullptr;
	delete[] helpIns;
	helpIns = nullptr;
	countHelpValues = 0;
	countLayers = 0;
	countIns = 0;
	data = nullptr;
}
// Связывает сеть с данными обучения. Сеть работает с данными, но не может их менять и удалять. Генерирует исключение, если:
// - countIns != _data.getCountIns();
// - topology[countLayers-1] != _data.getCountOuts().
void PWeb::bindData(const PData *_data)
{
	if (countIns != _data->getCountIns()) throw(std::exception("countIns != _data.getCountIns()"));
	if (topology[countLayers-1] != _data->getCountOuts()) throw(std::exception("topology[countLayers-1] != _data.getCountOuts()"));
	data = _data;
}
// Запускает обучение, используя вспомогательный буфер buffer. 
// _Time - максимальное время работы в миллисекундах;
// _maxIterations - максимальнок количество итераций;
// _sumSigma - минимальная ошибка, при достижении которой обучение заканчивается;
// _mForce - сила мутации;
// _nForce - сила шума при создании первичной популяции;
// _sForce - сила супермутации.
// Ничего не делает, если:
// - _countParents < 100;
// - _countIterations < 1;
// - _countParents < numCPU;
// - _numCPU > MAXIMUM_WAIT_OBJECTS.
// Обучение идет до тех пор, пока: 
// - не нажмут Escape, или не пройдет заданное кол-во итераций, или sumSigma сети не станет <= _sumSigma, 
// или пока не пройдет _maxTime миллисекунд.
// При _print == true после каждой итерации печатает статистику обучения.
// Не генерирует исключений.
void PWeb::learn(const size_t _countParents, const size_t _maxIterations, const double _sumSigma, const size_t _time, 
	const double _nForce, const double _mForce, const double _sForce, const bool _print, const size_t _numCPU) _NOEXCEPT
{
	SYSTEM_INFO systemInfo;
	GetSystemInfo(&systemInfo);
	size_t numCPU;
	if ((_countParents < 100) || (_maxIterations < 1) || (_numCPU > _countParents) || (_numCPU > MAXIMUM_WAIT_OBJECTS)) return;
	if (_numCPU == 0)
	{
		numCPU = systemInfo.dwNumberOfProcessors;
	}
	else {
		numCPU = _numCPU;
	}
	size_t time;
	if (_time == 0)
	{
		time = -1;
	}
	else {
		time = _time;
	}
	HANDLE *hStartSemaphores = new HANDLE[numCPU];
	HANDLE *hFinishSemaphores = new HANDLE[numCPU];
	atomic<bool> threadEnd{ false };
	// Готовим буфер для работы генетического алгоритма. 
	PWeb **buffer = new PWeb*[_countParents * _countParents];
	for (size_t w = 0; w < _countParents * _countParents; ++w)
	{
		buffer[w] = new PWeb(*this);
	}
	// Готовим массив с параметрами для потоков.
	ThreadParam *params = new ThreadParam[numCPU];
	HANDLE *hThreads = new HANDLE[numCPU];

	// Всех родителей, кроме _buffer[0] и _buffer[_countParents - 1], заполняем шумом.
	for (size_t w = 1; w < _countParents - 2; ++w)
	{
		buffer[w]->fillNoise(_nForce);
	}
	// Готовим параметры для рабочих потоков. И запускаем их.
	for (size_t t = 0; t < numCPU; ++t)
	{
		hFinishSemaphores[t] = CreateSemaphore(nullptr, 0, 1, nullptr);
		hStartSemaphores[t] = CreateSemaphore(nullptr, 0, 1, nullptr);
		params[t].numCPU = numCPU;
		params[t].ID = t;
		params[t].buffer = buffer;
		params[t].countParents = _countParents;
		params[t].mForce = _mForce;
		params[t].hStartSemaphore = hStartSemaphores[t];
		params[t].hFinishSemaphore = hFinishSemaphores[t];
		params[t].threadEnd = &threadEnd;
		hThreads[t] = (HANDLE)_beginthreadex(nullptr, 0, crossAndTestThread, &params[t], 0, 0);
	}
	// Цикл скрещивания и отбора.
	size_t iter = 0;
	size_t t0 = clock();
	do
	{
		// Последний предок становится жертвой супермутации.
		buffer[_countParents - 1]->fillNoise(_sForce);
		// Во множестве потоков скрещиваем и пропускаем через всех наследников уроки, определяя их качество.
		for (size_t t = 0; t < numCPU; ++t)
		{
			ReleaseSemaphore(hStartSemaphores[t], 1, 0);
		}
		WaitForMultipleObjects(numCPU, hFinishSemaphores, true, INFINITE);

		// Отсортируем потомков по величинам ошибок при помощи вспомогательного массива.
		vector<size_t> helpBuffer(_countParents * (_countParents - 1));// Можно оптимизировать ценой читабельности.
		for (size_t i = 0; i < helpBuffer.size(); ++i)
		{
			helpBuffer[i] = i + _countParents;
		}

		// Сортировка занимает < 0.1 % общего времени. Смысла в многопоточной сортировке нет.
		// Сортируем по sumSigma.
		sort(helpBuffer.begin(), helpBuffer.end(), [&buffer](const size_t &a, const size_t &b)->bool
		{
			return (buffer[a]->getSumSigma() < buffer[b]->getSumSigma());
		});
		// Отбираем лучших наследников.
		for (size_t k = 0; k < _countParents; ++k)
		{
			std::swap(buffer[k], buffer[helpBuffer[k]]);
		}

		if (_print == true)
		{
			cout << "Time elapsed: " << (clock() - t0) / 1000.f << " seconds" << endl;
			cout << "Iteration: " << iter << endl;
			cout << "sumSigma: " << buffer[0]->getSumSigma() << endl << endl;
		}
		// Повторяем скрещивание.
		++iter;
	} while ( !GetAsyncKeyState(VK_ESCAPE) && (iter < _maxIterations) && (buffer[0]->sumSigma >= _sumSigma) && (time >= (clock() - t0)) );
	*this = *buffer[0];
	// Освобождаем ресурсы.
	threadEnd.store(true);
	for (size_t t = 0; t < numCPU; ++t)
	{
		ReleaseSemaphore(hStartSemaphores[t], 1, 0);
	}
	WaitForMultipleObjects(numCPU, hFinishSemaphores, true, INFINITE);
	for (size_t t = 0; t < numCPU; ++t)
	{
		CloseHandle(hThreads[t]);
	}
	delete[] hThreads;
	hThreads = nullptr;
	for (size_t t = 0; t < numCPU; ++t)
	{
		CloseHandle(hStartSemaphores[t]);
		CloseHandle(hFinishSemaphores[t]);
	}
	delete[] hStartSemaphores;
	hStartSemaphores = nullptr;
	delete[] hFinishSemaphores;
	hFinishSemaphores = nullptr;
	delete[] params;
	params = nullptr;
	for (size_t w = 0; w < _countParents * _countParents; ++w)
	{
		delete buffer[w];
		buffer[w] = nullptr;
	}
	delete[] buffer;
	buffer = nullptr;
}
// Задаем сети вопрос. _insOuts - это массив, первые this->countIns значений - являются входными(вопросными).
// Сеть кладет ответы с выходов в последние this->countOuts ячеек. Размерность не проверяется. Не генерирует исключений.
void PWeb::ask(double *_insOuts) _NOEXCEPT
{
	// Пропускаем входные сигналы урока через нейроны первого слоя.
	for (size_t n = 0; n < topology[0]; ++n)
	{
		neurons[0][n]->genOut(_insOuts, helpOuts[n], f);
	}
	// Идем по слоям и пропускаем через нейроны каждого слоя выходные сигналы нейронов предыдущего слоя.
	for (size_t l = 1; l < countLayers; ++l)
	{
		std::swap(helpIns, helpOuts);
		for (size_t n = 0; n < topology[l]; ++n)
		{
			neurons[l][n]->genOut(helpIns, helpOuts[n], f);
		}
	}
	for (size_t n = 0; n < topology[countLayers - 1]; ++n)
	{
		_insOuts[n + countIns] = helpOuts[n];
	}
}
// Сохраняем сеть в бинарный файл. Выходные сигналы нейронов, функция активации нейронов и данные обучения не сохраняются. 
// Генерирует исключение, если:
// - _fileName = nullptr.
void PWeb::saveToBinFile(const char *_fileName)
{
	fstream f;
	if (_fileName == nullptr) throw(std::exception("_fileName == nullptr"));
	f.open(_fileName, ios_base::out | ios_base::binary);
	f.write((const char*)&countIns,sizeof(size_t));
	f.write((const char*)&countLayers, sizeof(size_t));
	f.write((const char*)topology, sizeof(size_t) * countLayers);
	f.write((const char*)&countHelpValues, sizeof(size_t));
	for (size_t l = 0; l < countLayers; ++l)
	{
		for (size_t n = 0; n < topology[l]; ++n)
		{
			f.write( (const char*)&neurons[l][n]->getCountWeights(), sizeof(size_t));
			f.write((const char*)&neurons[l][n][0][0], sizeof(double) * neurons[l][n]->getCountWeights());
		}
	}
	f.write((const char*)&sumSigma, sizeof(double));
	f.close();
}
// Загружаем сеть из бинарного файла. Выходные сигналы нейронов, helpIns, helpOuts, функции активации нейронов и данные обучения не загружаются.
// Генерирует исключение, если:
// _fileName = nullptr;
// - файл не найден.
void PWeb::loadFromBinFile(const char *_fileName)// Надо бы сделать конструктор загрузки из бинарного файла.
{
	fstream f;
	if (_fileName == nullptr) throw(std::exception("_fileName == nullptr"));
	f.open(_fileName, ios_base::in | ios_base::binary);
	if (f.is_open() == false) throw(std::exception("No such file"));
	// Меняем нашу топологию только в том случае, если в этом есть необходимость. Но меняем полностью.
	bool equal = true;
	size_t helpCountIns{ 0 };
	size_t helpCountLayers{ 0 };
	size_t *helpTopology{ 0 };
	f.read((char*)&helpCountIns, sizeof(size_t));
	f.read((char*)&helpCountLayers, sizeof(size_t));
	helpTopology = new size_t[helpCountLayers]; 
	f.read((char*)helpTopology, sizeof(size_t)*helpCountLayers);
	size_t helpCountHelpValues{ 0 };
	f.read((char*)&helpCountHelpValues, sizeof(size_t));
	if (helpCountIns != countIns)
	{
		equal = false;
		goto next;
	}
	for (size_t l = 0; l < countLayers; ++l)
	{
		if (helpTopology[l] != topology[l])
		{
			equal = false;
			goto next;
		}
	}
	next:
	if (equal == true)
	{
		for (size_t l = 0; l < countLayers; ++l)
		{
			for (size_t n = 0; n < topology[l]; ++n)
			{
				f.read((char*)&neurons[l][n]->getCountWeights(), sizeof(size_t));
				f.read((char*)&neurons[l][n][0][0], sizeof(double) * neurons[l][n]->getCountWeights());
			}
		}
	} else {
		this->~PWeb();
		countIns = helpCountIns;
		countLayers = helpCountLayers;
		std::swap(topology, helpTopology);
		countHelpValues = helpCountHelpValues;
		helpIns = new double[countHelpValues];
		helpOuts = new double[countHelpValues];
		neurons = new PNeuron**[countLayers];
		for (size_t l = 0; l < countLayers; ++l)
		{
			neurons[l] = new PNeuron*[topology[l]];
			for (size_t n = 0; n < topology[l]; ++n)
			{
				size_t helpCountWeights{ 0 };
				f.read((char*)&helpCountWeights, sizeof(size_t));
				neurons[l][n] = new PNeuron(helpCountWeights);
				f.read((char*)&neurons[l][n][0][0], sizeof(double) * helpCountWeights);
			}
		}
	}
	f.read((char*)&sumSigma, sizeof(double));
	f.close();
	delete[] helpTopology;
}

// Тестируем сеть. Пропускаем через нее все уроки, чтобы определить ошибки. Ничего не делает, если data == nullptr.
// Не генерирует исключений.
void PWeb::test() _NOEXCEPT
{
	if (data == nullptr) return;
	sumSigma = 0.0;
	size_t countValuesInLesson = data->getCountIns() + data->getCountOuts();
	for (size_t e = 0; e < data->getCountLessons(); ++e)
	{
		// Пропускаем входные сигналы урока через нейроны первого слоя.
		memcpy(helpIns, &data->getData()[e * countValuesInLesson], sizeof(double)* data->getCountIns());
		for (size_t n = 0; n < topology[0]; ++n)
		{
			neurons[0][n]->genOut(helpIns, helpOuts[n], f);
		}
		// Идем по слоям и пропускаем через нейроны каждого слоя выходные сигналы нейронов предыдущего слоя.
		for (size_t l = 1; l < countLayers; ++l)
		{
			std::swap(helpIns, helpOuts);
			for (size_t n = 0; n < topology[l]; ++n)
			{
				neurons[l][n]->genOut(helpIns, helpOuts[n], f);
			}
		}
		// Т.к. вычисляется абсолютная ошибка - в первую очередь аппроксимация будет происходить по наибольшим выходным сигналам.
		// Но, при sumSigma -> 0, отклонения слабых выходных сигналов от требуемых будут -> 0.
		for (size_t o = 0; o < topology[countLayers - 1]; ++o)
		{
			sumSigma += abs(helpOuts[o] - data->getData()[e * countValuesInLesson + countIns + o]);
		}
	}
}
// Возвращает ссылку на sumSigma. Не генерирует исключений.
const double &PWeb::getSumSigma() const _NOEXCEPT
{
	return sumSigma;
}
// Заполняет веса нейронаов сети случайными величинами с равномерным распределением из интервала [-1.0 * _nForce; +1.0 * _nForce). 
// Не генерирует исключений.
void PWeb::fillNoise(const double _nForce) _NOEXCEPT
{
	for (size_t l = 0; l < countLayers; ++l)
	{
		for (size_t n = 0; n < topology[l]; ++n)
		{
			for (size_t w = 0; w < neurons[l][n]->getCountWeights(); ++w)
			{
				neurons[l][n][0][w] = globalWDistribution(globalRandomEngine) * _nForce;
			}
		}
	}
}
// Конструктор по умолчанию. Не генерирует исключений.
PData::PData() _NOEXCEPT : countIns{ 0 }, countOuts{ 0 }, countLessons{ 0 }, data{ nullptr }
{
}
// Конструктор копирования. Не генерирует исключений.
PData::PData(const PData &_object) _NOEXCEPT : countIns{ _object.countIns }, countOuts{ _object.countOuts }, countLessons{ _object.countLessons },
data{ nullptr }
{
	if (countLessons != 0)
	{
		data = new double[countLessons * (countIns + countOuts)];
		memcpy(data, _object.data, sizeof(double)*(countLessons*(countIns+countOuts)));
	}
}
// Конструктор перемещения. Не генерирует исключений.
PData::PData(PData &&_object) _NOEXCEPT : countIns{ _object.countIns }, countOuts{ _object.countOuts }, countLessons{ _object.countLessons },
data{ _object.data }
{
	_object.countIns = 0;
	_object.countOuts = 0;
	_object.countLessons = 0;
	_object.data = nullptr;
}
// Копирующий оператор присваивания, полностью оптимизирован, перестраивает data только при необходимости, 
// но перестраивает полностью. Не генерирует исключений.
PData &PData::operator=(const PData &_object) _NOEXCEPT
{
	if (this != &_object)// Подумать, что будет, если countRecords то же, но countIns/Outs различаются!
	{
		if (countLessons == _object.countLessons)
		{
			countIns = _object.countIns;
			countOuts = _object.countOuts;
			memcpy(data, _object.data, sizeof(double)*(countLessons*(countIns+countOuts)));
		}
		else {
			countIns = _object.countIns;
			countOuts = _object.countOuts;
			if (countLessons != 0) delete[] data;
			countLessons = _object.countLessons;
			data = new double[countLessons*(countIns+countOuts)];
			memcpy(data, _object.data, sizeof(double)*(countLessons*(countIns+countOuts)));
		}
	}
	return *this;
}
// Перемещающий оператор присваивания. Не генерирует исключений.
PData &PData::operator=(PData &&_object) _NOEXCEPT
{
	if (this != &_object)
	{
		delete[] data;
		countIns = _object.countIns;
		countOuts = _object.countOuts;
		countLessons = _object.countLessons;
		data = _object.data;
		_object.countIns = 0;
		_object.countOuts = 0;
		_object.countLessons = 0;
		_object.data = nullptr;
	}
	return *this;
}
// Деструктор с полным занулением. Не генерирует исключений.
PData::~PData() _NOEXCEPT
{
	countIns = 0;
	countOuts = 0;
	countLessons = 0;
	delete[] data;
	data = nullptr;
}
// Полное очищение. Не генерирует исключений.
void PData::clear() _NOEXCEPT
{
	this->~PData();
}
// Загрузка данных из текстового файла. Возвращает число загруженных уроков. Если данные из файла удается загрузить, то все старые
// данные стираются. Генерирует исключение, если:
// - файл не найден;
// - _countIns < 1;
// - _countOuts < 1;
// - число записей в файле не кратно (_countIns + _countOuts).
size_t PData::loadFromTxtFile(const size_t _countIns, const size_t _countOuts, const const char *_fileName)
{
	fstream f;
	f.open(_fileName, ios_base::in);
	if (f.is_open() == false)
	{
		throw(std::exception("No such file"));
	}
	if (_countIns < 1)
	{
		throw(std::exception("_countIns < 1"));
	}
	if (_countOuts < 1) 
	{
		throw(std::exception("_countOuts < 1"));
	}
	double helpValue;
	size_t countValues = 0;
	while (!f.eof())
	{
		helpValue = numeric_limits<double>::signaling_NaN();
		f >> helpValue;
		if (helpValue != numeric_limits<double>::signaling_NaN())
		{
			++countValues;
		}
	}
	
	if ((countValues % (_countIns + _countOuts)) != 0)
	{
		f.close();
		throw(std::exception("countValues % (_countIns + _countOuts) != 0"));
	}
	delete[] data;
	countIns = _countIns;
	countOuts = _countOuts;
	countLessons = countValues / (_countIns + _countOuts);
	data = new double[countValues];
	f.seekg(0, ios::beg);
	size_t d = 0;
	while (!f.eof())
	{
		f >> helpValue;
		data[d] = helpValue;
		++d;
	}
	f.close();
	return countLessons;
}
// Возвращает указатель на data. Не генерирует исключений.
const double *PData::getData() const _NOEXCEPT
{
	return data;
}
// Возвращает ссылку countLessons. Не генерирует исключений.
const size_t &PData::getCountLessons() const _NOEXCEPT
{
	return countLessons;
}
// Возвращает ссылку на константный countIns. Не генерирует исключений.
const double &PData::getCountIns() const _NOEXCEPT
{
	return countIns;
}
// Возвращает ссылку на константный countOuts. Не генерирует исключений.
const double &PData::getCountOuts() const _NOEXCEPT
{
	return countOuts;
}
// Выводит в консоль все данные.
void PData::printData() const _NOEXCEPT
{
	cout << "Lessons count: " << countLessons << endl;
	for (size_t l = 0; l < countLessons; ++l)
	{
		cout << "Lesson_" << l << " Ins: ";
		for (size_t i = 0; i < countIns; ++i)
		{
			cout << data[l*(countIns + countOuts) + i] << " ";
		}
		cout << "Outs: ";
		for (size_t o = 0; o < countOuts; ++o)
		{
			cout << data[l*(countIns + countOuts) + countIns + o] << " ";
		}
		cout << endl;
	}
}
// Функция для потока, в котором будет происходить скрещивание и тестирование потомков. Не генерирует исключений.
unsigned int _stdcall PWeb::crossAndTestThread(void *_param) _NOEXCEPT
{
	default_random_engine threadRandomEngine( clock() );
	uniform_real_distribution<double> threadWDistribution(globalWDistribution);
	uniform_int_distribution<size_t> threadProbobilityCrossDistribution(globalProbobilityCrossDistribution);
	uniform_int_distribution<size_t> threadProbobilityMutationDistribution(globalProbobilityMutationDistribution);
	const size_t numCPU = ((ThreadParam*)_param)->numCPU;
	const size_t ID = ((ThreadParam*)_param)->ID;
	PWeb **buffer = ((ThreadParam*)_param)->buffer;
	const size_t countParents = ((ThreadParam*)_param)->countParents;
	const double mForce = ((ThreadParam*)_param)->mForce;
	const HANDLE hStartSemaphore = ((ThreadParam*)_param)->hStartSemaphore;
	const HANDLE hFinishSemaphore = ((ThreadParam*)_param)->hFinishSemaphore;
	const atomic<bool> *threadEnd = ((ThreadParam*)_param)->threadEnd;
	const size_t countLayers = buffer[0]->countLayers;
	size_t *topology = new size_t[countLayers];
	memcpy(topology, buffer[0]->topology, sizeof(size_t) * countLayers);
	size_t startParent;
	size_t finishParent;
	// Адаптивная подстройка под различные значения countParents и numCPU.
	size_t startOurChilds{ countParents };
	size_t helpValue;
	for (size_t i = 0; i < ID; ++i)
	{	
		helpValue = 0;
		for (size_t j = i; j < countParents; j += numCPU)
		{
			++helpValue;
		}
		startOurChilds += helpValue * (countParents - 1);
	}
	// Ожидаем команды.
	A:
	ptrdiff_t wResult = WaitForSingleObject(hStartSemaphore, INFINITE);
	switch (wResult)
	{
		case(WAIT_OBJECT_0) :
		{
			if (threadEnd->load() == false)
			{
				// Совершаем вычисления и возвращаемся к ожиданию.
				size_t localBias = 0;
				double mutation;
				size_t stepsToNextMutation = 1;
				for (size_t p1 = ID; p1 < countParents; p1 += numCPU)
				{
					// Т.к. множество потоков работают с предками, изменять их в этот момент нельзя.
					// Мутации подвергаются связи наследников. Непосредственно в момент скрещивания.
					for (size_t p2 = 0; p2 < countParents; ++p2)// Поправить хвостик! Поток, которому достался хвостик, выведет больше.
					{
						if (p1 != p2)
						{
							// Скрещиваем.
							for (size_t l = 0; l < countLayers; ++l)
							{
								for (size_t n = 0; n < topology[l]; ++n)
								{
									for (size_t w = 0; w < buffer[0]->neurons[l][n]->getCountWeights(); ++w)
									{
										// Мутация.
										mutation = 0.0;
										// Оптимизация (т.к. мутация - это довольно редкое явление).
										if (--stepsToNextMutation == 0)
										{
											stepsToNextMutation = threadProbobilityMutationDistribution(threadRandomEngine);
											mutation = threadWDistribution(threadRandomEngine) * mForce;
										}
										// Скрещивание.
										if (threadProbobilityCrossDistribution(threadRandomEngine) <= 50)
										{
											buffer[startOurChilds + localBias]->neurons[l][n][0][w] =
													buffer[p1]->neurons[l][n][0][w] + mutation;
										}
										else {
											buffer[startOurChilds + localBias]->neurons[l][n][0][w] =
													buffer[p2]->neurons[l][n][0][w] + mutation;
										}
									}
								}
							}
							// Тестируем.
							buffer[startOurChilds + localBias]->test();
							++localBias;
						}
					}
				}
				ReleaseSemaphore(hFinishSemaphore, 1, 0);
				goto A;
			}
			else {
				goto B;
			}
		}
		break;
	}
	B:
	delete[] topology;
	topology = nullptr;
	ReleaseSemaphore(hFinishSemaphore, 1, 0);
	return 0;
}
// Глобальный генератор случайных чисел.
default_random_engine GP::globalRandomEngine( clock() );
// Равномерное распределение для генерации случайной величины из интервала [-1.0; +1.0).
uniform_real_distribution<double> GP::globalWDistribution(-1, +1);
// Равномерное распределение для скрещивания.
uniform_int_distribution<size_t> GP::globalProbobilityCrossDistribution(1, 100);
// Равномерное распределение для мутации.
uniform_int_distribution<size_t> GP::globalProbobilityMutationDistribution(1, 200);
