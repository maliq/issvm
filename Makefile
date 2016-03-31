# options are "debug", "optimize" and "profile"
BUILD := debug

AR   := ar
MAKE := make
CC   := gcc
CXX  := g++
LD   := g++


#====  files  =================================================================

TARGET_SOURCES := \
	issvm_evaluate.cpp \
	issvm_initialize.cpp \
	issvm_optimize.cpp \
	issvm_recalculate.cpp

HEADERS := \
	array_sum.hpp \
	data.hpp \
	fast_exp.hpp \
	random_distribution_discrete_uniform.hpp \
	random_distribution.hpp \
	random_distribution_private_ziggurat.hpp \
	random_distribution_standard_exponential.hpp \
	random_distribution_standard_gamma.hpp \
	random_distribution_standard_gaussian.hpp \
	random_distribution_standard_uniform.hpp \
	random_generator.hpp \
	random_generator_lagged_fibonacci_4.hpp \
	random_generator_linear_congruential.hpp \
	random_generator_private_helpers.hpp \
	random_generator_private_lagged_fibonacci_4_helper.hpp \
	random_generator_private_linear_congruential_helper.hpp \
	random.hpp \
	svm.hpp \
	svm_kernel_base.hpp \
	svm_kernel_construct.hpp \
	svm_kernel.hpp \
	svm_kernel_private_cache.hpp \
	svm_kernel_private_data.hpp \
	svm_kernel_simple.hpp \
	svm_kernel_traits_gaussian.hpp \
	svm_kernel_traits.hpp \
	svm_kernel_traits_linear.hpp \
	svm_kernel_vector_data.hpp \
	svm_optimizer_base.hpp \
	svm_optimizer_classification_biased_perceptron.hpp \
	svm_optimizer_classification_biased_sbp.hpp \
	svm_optimizer_classification_biased_smo.hpp \
	svm_optimizer_classification_biased_sparsifier.hpp \
	svm_optimizer_classification_private_find_water_level.hpp \
	svm_optimizer_classification_unbiased_perceptron.hpp \
	svm_optimizer_classification_unbiased_sbp.hpp \
	svm_optimizer_classification_unbiased_smo.hpp \
	svm_optimizer_classification_unbiased_sparsifier.hpp \
	svm_optimizer.hpp \
	svm_serialization.hpp \
	vector.hpp

SOURCES := \
	svm_kernel_base.cpp \
	svm_kernel_private_cache.cpp \
	svm_optimizer_base.cpp \
	svm_optimizer_classification_biased_perceptron.cpp \
	svm_optimizer_classification_biased_sbp.cpp \
	svm_optimizer_classification_biased_smo.cpp \
	svm_optimizer_classification_biased_sparsifier.cpp \
	svm_optimizer_classification_private_find_water_level.cpp \
	svm_optimizer_classification_unbiased_perceptron.cpp \
	svm_optimizer_classification_unbiased_sbp.cpp \
	svm_optimizer_classification_unbiased_smo.cpp \
	svm_optimizer_classification_unbiased_sparsifier.cpp

LIBRARIES := \
	stdc++ \
	m

BOOST_LIBRARIES := \
	iostreams \
	program_options \
	regex \
	serialization


#====  platform-dependent  ====================================================

DEFINE_FLAGS :=

COMPILER_SHARED_FLAGS   := -D_GNU_SOURCE -Wall -fPIC -msse2 -mfpmath=sse -fopenmp
COMPILER_DEBUG_FLAGS    := -g3
COMPILER_OPTIMIZE_FLAGS := -O2 -funroll-loops -fomit-frame-pointer -DNDEBUG
COMPILER_PROFILE_FLAGS  := -g3 -fno-inline $(COMPILER_OPTIMIZE_FLAGS)

LINKER_FLAGS := \
	-fopenmp \
	${patsubst %,-l%,$(LIBRARIES)} \
	${patsubst %,-lboost_%,$(BOOST_LIBRARIES)}


#====  platform-independent  ==================================================

ifeq ($(BUILD),debug)
COMPILER_FLAGS := \
	$(DEFINE_FLAGS) \
	$(COMPILER_SHARED_FLAGS) \
	$(COMPILER_DEBUG_FLAGS)
endif

ifeq ($(BUILD),optimize)
COMPILER_FLAGS := \
	$(DEFINE_FLAGS) \
	$(COMPILER_SHARED_FLAGS) \
	$(COMPILER_OPTIMIZE_FLAGS)
endif

ifeq ($(BUILD),profile)
COMPILER_FLAGS := \
	$(DEFINE_FLAGS) \
	$(COMPILER_SHARED_FLAGS) \
	$(COMPILER_PROFILE_FLAGS)
endif

LDFLAGS := ${strip $(LINKER_FLAGS)}
CFLAGS := ${strip -std=gnu99 $(COMPILER_FLAGS) -I /Users/ricardo/boost/boost_1_58_0}
CXXFLAGS := ${strip -std=gnu++0x $(COMPILER_FLAGS) -I /Users/ricardo/boost/boost_1_58_0}


#====  derived variables  =====================================================

OBJECTS := \
	${patsubst %.cpp,%.o,${filter %.cpp,$(SOURCES)}} \
	${patsubst %.c,%.o,${filter %.c,$(SOURCES)}}
OBJECTS := ${strip $(OBJECTS)}

TARGET_OBJECTS := \
	${patsubst %.cpp,%.o,${filter %.cpp,$(TARGET_SOURCES)}} \
	${patsubst %.c,%.o,${filter %.c,$(TARGET_SOURCES)}}
TARGET_OBJECTS := ${strip $(TARGET_OBJECTS)}

TARGETS := \
	${patsubst %.cpp,%,${filter %.cpp,$(TARGET_SOURCES)}} \
	${patsubst %.c,%,${filter %.c,$(TARGET_SOURCES)}}
TARGETS := ${strip $(TARGETS)}


#====  compilation rules  =====================================================

# disable default rules
.SUFFIXES:

.PHONY : all
all: $(OBJECTS) $(TARGET_OBJECTS) $(TARGETS)

%.o : %.c $(HEADERS)
	@echo "----  Building \"$@\" from \"$<\"  ----"
	$(CC) -c $< -o $@ $(CFLAGS)
	@echo

%.o : %.cpp $(HEADERS)
	@echo "----  Building \"$@\" from \"$<\"  ----"
	$(CXX) -c $< -o $@ $(CXXFLAGS)
	@echo

% : %.o $(OBJECTS)
	@echo "----  Building \"$@\" from \"$<\"  ----"
	$(LD) $< $(OBJECTS) -o $@ $(LDFLAGS)
	@echo

.PHONY : clean
clean :
	@echo "----  Cleaning  ----"
	rm -f $(OBJECTS) $(TARGET_OBJECTS) $(TARGETS)
	@echo
