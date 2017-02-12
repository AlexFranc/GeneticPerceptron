#ifndef GP_H
#define GP_H

#include <initializer_list>
using namespace std;

namespace GP
{
	// Предварительное объявление прототипов, для подсказок компилятору.
	class PData;
	class PNeuron;
	class PWeb;

	// Прототипы.

	// Параметр для потока тестирования.

	struct ThreadParam
	{
		size_t numCPU;
		size_t ID;
		PWeb **buffer;
		size_t countParents;
		double mForce;
		HANDLE hStartSemaphore;
		HANDLE hFinishSemaphore;
		atomic<bool> *threadEnd;
	};
	
	// PData потоко-небезопасен. Пока объект PWeb обучается, используя объект PData, запрещается как-либо изменять содержимое используемого 
	// объекта PData.

	class PData
	{
	public:
		PData() _NOEXCEPT;
		PData(const PData &_object) _NOEXCEPT;
		PData(PData &&_object) _NOEXCEPT;
		PData &operator=(const PData &_object) _NOEXCEPT;  
		PData &operator=(PData &&_object) _NOEXCEPT;
		~PData() _NOEXCEPT;
		void clear() _NOEXCEPT;
		size_t loadFromTxtFile(const size_t _countIns, const size_t _countOuts, const const char *_fileName);
		const double *getData() const _NOEXCEPT;
		const size_t &getCountLessons() const _NOEXCEPT;
		const double &getCountIns() const _NOEXCEPT;
		const double &getCountOuts() const _NOEXCEPT;
		void printData() const _NOEXCEPT;
	protected:
	private:
		size_t countIns;
		size_t countOuts;
		size_t countLessons;
		double *data;
	};
	class PWeb
	{
	public:
		PWeb() = delete;
		PWeb(const size_t _countIns, initializer_list<size_t> _topology, double(*_f)(double) = nullptr);
		PWeb(const PWeb &_object) _NOEXCEPT;
		PWeb(PWeb &&_object) _NOEXCEPT;
		PWeb &operator=(const PWeb &_object) _NOEXCEPT;
		PWeb &operator=(PWeb &&_object) _NOEXCEPT;
		~PWeb() _NOEXCEPT;
		void bindData(const PData *_data);
		void learn(const size_t _countParents = 100, const size_t _maxIterations = 1000, const double _sumSigma = 0.0, 
			const size_t _time = 0, const double _nForce = 0.1, const double _mForce = 1.0, const double _sForce = 1E9, 
			const bool _print = true, const size_t _numCPU = 0) _NOEXCEPT;
		void test() _NOEXCEPT;
		const double &getSumSigma() const _NOEXCEPT;
		void fillNoise(const double _nForce = 1.0) _NOEXCEPT;
		void ask(double *_insOuts) _NOEXCEPT;
		void saveToBinFile(const char *_fileName);
		void loadFromBinFile(const char *_fileName);
	protected:
	private:
		size_t countIns;
		size_t countLayers;
		size_t *topology;
		size_t countHelpValues;
		size_t countWeights;
		double *helpIns;
		double *helpOuts;
		double *weights;
		double sumSigma;
		double (*f)(double);
		const PData *data;
		
		static unsigned int _stdcall crossAndTestThread(void*) _NOEXCEPT;
	};

	extern default_random_engine globalRandomEngine;
	extern uniform_real_distribution<double> globalValueDistribution;
};
#endif
