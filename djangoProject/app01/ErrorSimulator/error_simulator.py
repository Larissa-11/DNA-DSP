import copy
import random
from scipy.stats import skewnorm
from tqdm import tqdm

from app01.ErrorSimulator import strand_error_simulation
from app01.ErrorSimulator.custom_random_variable import CustomRvContinuous
from django.views import View


def parse_rate(rate_str) -> float:
    if isinstance(rate_str, float):
        return rate_str
    index = rate_str.find('E')
    if index == -1:
        return float(rate_str)
    else:
        num = float(rate_str[:index])
        exp = 10 ** float(rate_str[index + 1:])
        return num * exp


def parse_rates_dictionary(rates_dict):
    for key, value in rates_dict.items():
        if isinstance(value, dict):
            for error, rate in value.items():
                rates_dict[key][error] = parse_rate(rate)
        else:
            rates_dict[key] = parse_rate(value)


def mess_output_strands(evyat_content):
    lines = evyat_content.split('\n')
    shuffled_lines = []

    iterator = iter(lines)

    next(iterator)
    next(iterator)

    for line in iterator:
        if line == '':
            # Skip current line (which is an empty string) and the next item (another empty string) + 2 items:
            # Here, we need to get the next three items from the iterator.
            try:
                next(iterator)  # Skip the first ''
                next(iterator)  # Skip the second empty line
                next(iterator)  # Skip the origin strand
                # The next iteration will point to the next output strand
            except StopIteration:
                break  # 如果没有更多元素，退出循环
        else:
            shuffled_lines.append(line)
    random.shuffle(shuffled_lines)
    return '\n'.join(shuffled_lines)


class Simulator(View):
    def __init__(self, total_error_rates, base_error_rates, min_copies, max_copies, input_file, is_stutter_method=False,
                 distribution_info=None):
        self.total_error_rates = copy.deepcopy(total_error_rates)
        self.base_error_rates = copy.deepcopy(base_error_rates)
        parse_rates_dictionary(self.total_error_rates)
        parse_rates_dictionary(self.base_error_rates)
        self.long_deletion_length_rates = {2: 2.8 * (10 ** (-4)),
                                           3: 7.75 * (10 ** (-5)),
                                           4: 3.25 * (10 ** (-5)),
                                           5: 10 ** (-6),
                                           6: 5.5 * (10 ** (-8))}
        self.input_file = input_file
        self.is_stutter_method = is_stutter_method
        self.distribution_info = distribution_info
        self.random = None
        self.min_copies = min_copies
        self.max_copies = max_copies - 1

    def simulate_errors(self):
        input_f = self.input_file.read().decode('utf-8').splitlines()
        num_values = len(input_f)

        if self.distribution_info is None:
            skewness = 10
            self.random = skewnorm.rvs(a=skewness, loc=self.max_copies, size=num_values)
            self.random = self.random - min(self.random)
            self.random = self.random / max(self.random)
            self.random = self.random * self.max_copies
            self.random = self.random + self.min_copies
            self.random = [round(x) for x in self.random]
        elif self.distribution_info['type'] == 'continuous':
            raw_samples = CustomRvContinuous.multithreaded_rvs(size=num_values,
                                                               pdf_str=self.distribution_info['value'],
                                                               min_value=self.distribution_info['min'],
                                                               max_value=self.distribution_info['max'])
            self.random = [round(x) for x in raw_samples]
        elif self.distribution_info['type'] == 'vector':
            self.random = self.distribution_info['value']
            self.min_copies = min(self.random)
            self.max_copies = max(self.random)

        evyat_content = []
        for i, line in enumerate(tqdm(input_f, desc='Simulation', ncols=100, total=num_values)):
            original_strand = line.rstrip()
            evyat_content.append(original_strand + '\n' + '*****************************')

            num_copies = self.min_copies
            if i < len(self.random):
                num_copies = self.random[i]
            else:
                num_copies = random.randint(self.min_copies, self.max_copies + 1)

            for j in range(num_copies):
                output_strand = copy.deepcopy(original_strand)
                strand_error_simulator = strand_error_simulation.StrandErrorSimulation(self.total_error_rates,
                                                                                       self.base_error_rates,
                                                                                       self.long_deletion_length_rates,
                                                                                       output_strand)
                if self.is_stutter_method:
                    output_strand = strand_error_simulator.simulate_stutter_errors_on_strand()
                else:
                    output_strand = strand_error_simulator.simulate_errors_on_strand()

                evyat_content.append(output_strand)

            evyat_content.append('\n')

        # 过滤掉空行
        # evyat_content = [line for line in evyat_content if line.strip()]

        evyat_str = '\n'.join(evyat_content)
        shuffled_str = mess_output_strands(evyat_str)

        return evyat_str, shuffled_str

    # def post(self, request):
    #     total_error_rates = request.POST.get('total_error_rates')
    #     base_error_rates = request.POST.get('base_error_rates')
    #     min_copies = int(request.POST.get('min_copies'))
    #     max_copies = int(request.POST.get('max_copies'))
    #     is_stutter_method = request.POST.get('is_stutter_method', 'false').lower() == 'true'
    #     distribution_info = request.POST.get('distribution_info', None)
    #
    #     if distribution_info:
    #         distribution_info = json.loads(distribution_info)
    #
    #     input_file = request.FILES['input_file']
    #
    #     sim = Simulator(total_error_rates, base_error_rates, min_copies, max_copies, input_file, is_stutter_method,
    #                     distribution_info)
    #     evyat_str, shuffled_str = sim.simulate_errors()
    #
    #     return JsonResponse({'evyat': evyat_str, 'shuffled': shuffled_str})
