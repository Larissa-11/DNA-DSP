import random
import copy


class StrandErrorSimulation:
    def __init__(self, total_error_rates, base_error_rates, deletion_length_rates, strand):
        """
        :param total_error_rates: Dictionary of the total error rates used in the simulation.
            Example of a dictionary:
            {'d': 0.1, 'i': 0.2, 's': 0.1, 'ld': 0.6}
        :param base_error_rates: Dictionary of dictionaries for each base.
            Example:
            {   'A': {'s': 0.1, 'i': 0.2, 'pi': 0.1, 'd': 0.05, 'ld': 0.6},
                'T': {...},
                'C': {...},
                'G': {...}
            }
        :param deletion_length_rates: Dictionary of deletion length rates for higher lengths than 1 (start from 2):
            Example:
            {2: 0.2, 3: 0.1, 4: 0.3, 5: 0.05, 6: 0.001}
        :param strand: The strand to implement errors on.
        """
        self.total_error_rates = total_error_rates
        self.base_error_rates = base_error_rates
        self.deletion_length_rates = deletion_length_rates
        self.strand = strand
        self.index = 0
        # for testing only:
        self.err_type = None

    ''' Main Methods: '''

    def simulate_stutter_errors_on_strand(self) -> str:
        """
        Simulates stutter errors on the given strand and returns the target strand.
        Use only for stutter.
        :return:
        Modified strand after errors simulation
        """
        while self.index < len(self.strand):
            self.simulate_stutter_error_on_base()
            self.index += 1
        # after implementing the stutter errors, go over the strand again and inject substitutions according to rates:
        self.index = 0
        while self.index < len(self.strand):
            base = self.strand[self.index]
            base_substitution_rate = self.base_error_rates[base]['s']
            options = ['y', 'n']
            rates = [base_substitution_rate, 1 - base_substitution_rate]
            draw = random.choices(options, weights=rates, k=1)
            if draw[0] == 'y':
                self.strand = self.inject_substitution()
            self.index += 1
        return self.strand

    def simulate_errors_on_strand(self) -> str:
        """
        Simulates errors on the given strand and returns the target strand.
        Use for any method EXCEPT stutter.
        :return:
        Modified strand after errors simulation
        """
        while self.index < len(self.strand):
            self.simulate_error_on_base()
            self.index += 1
        return self.strand

    ''' Helper Methods: '''

    def simulate_stutter_error_on_base(self):
        """
        Simulates a stutter error on the current base (in the current index)
        Modifies the working strand stored in class.
        NO RETURN VALUE
        """
        # this method doesn't have long deletion
        base = self.strand[self.index]
        options = ['y', 'n']
        base_stutter_rate = self.base_error_rates[base]['i']
        rates = [base_stutter_rate, 1 - base_stutter_rate]
        draw = random.choices(options, weights=rates, k=1)
        is_stutter = (draw[0] == 'y')
        # for testing, count stutter times
        self.err_type = 0
        while is_stutter:
            self.strand = self.strand[:self.index] + base + self.strand[self.index:]
            self.index += 1
            # for testing, count stutter times
            self.err_type += 1
            # draw again before next iteration:
            draw = random.choices(options, weights=rates, k=1)
            is_stutter = (draw[0] == 'y')
        self.strand = self.inject_error('d')
        return

    def simulate_error_on_base(self):
        """
        Simulates any error EXCEPT stutter on the current base (in the current index).
        Modifies the working strand stored in class.
        """
        base = self.strand[self.index]
        # 1. summarize all error rates into total rate, and conclude the complementary non-error rate:
        # 检查序列的总错误率不能超过1
        total_error_rate = 0
        for value in self.total_error_rates.values():
            total_error_rate += value
        no_error_rate = 1 - total_error_rate
        assert (no_error_rate >= 0)

        # 2. draw whether there's error or not in the given rates:
        # Give higher rates to first third
        options = ['y', 'n']
        rates = [total_error_rate, no_error_rate]
        if self.index <= (1 / 3) * len(self.strand):
            rates[0] = rates[0] * 3 / 2
            rates[1] = 1 - rates[0]
        else:
            rates[0] = rates[0] * 3 / 4
            rates[1] = 1 - rates[0]
        draw = random.choices(options, weights=rates, k=1)

        # 3. check type of drawn result:
        # 3.1. If there's error:
        if draw[0] == 'y':
            # generate an error type:
            error_type = self.generate_error_type_for_base(base)
            self.err_type = error_type  # for testing only
            self.strand = self.inject_error(error_type)
        # 3.2. If there's no error - do nothing.
        else:
            self.err_type = 'n'  # for testing only

    def generate_error_type_for_base(self, base) -> str:
        """
        Generate an error from the error rates dictionary passed as arguments:
        Returns the error type generated for the base (string).
        :param base: The value of the base currently working on: 'A', 'T', 'C', 'G'.
        :return 'd': for deletion
        :return 'ld': for long deletion
        :return 'pi': for insertion (base on pre-insertion symbol)
        :return 's': for substitution
        """
        # create two lists of the dictionary - options list and rates list:
        options = []
        rates = []
        base_rates = copy.deepcopy(self.base_error_rates[base])

        # remove insertion rates (as they are not needed in this stage - they are used for inserted base generation)
        del base_rates['i']

        for key, value in base_rates.items():
            options.append(key)
            rates.append(value)

        #  Note: choices uses weights, and thus is equivalent to conditional probability!
        draw = random.choices(options, weights=rates, k=1)
        # return error type string:
        return draw[0]

    def inject_deletion(self, error_type) -> str:
        """
        Inject deletion to the given strand, starting from the base in `index` location of the strand.
        Returns a strand with the injected error.
        :param error_type: Type of deletion: 'd' for single base deletion or 'ld' for long deletion.
        :return: a strand with the injected deletion.
        """
        modified_strand = ""

        if error_type == 'd':
            # single base deletion:
            if self.index == len(self.strand) - 1:
                modified_strand = self.strand[:self.index]
            else:
                modified_strand = self.strand[:self.index] + self.strand[self.index + 1:]

        elif error_type == 'ld':
            # multiple base deletion:
            long_del_dict = copy.deepcopy(self.deletion_length_rates)

            # draw length:
            options = list(long_del_dict.keys())
            rates = list(long_del_dict.values())
            draw = random.choices(options, weights=rates, k=1)

            deletion_length = draw[0]
            if self.index + deletion_length > len(self.strand) - 1:
                modified_strand = self.strand[:self.index]
            else:
                modified_strand = self.strand[:self.index] + self.strand[self.index + deletion_length:]

        # keep index the same! The original base in that index (or further) was deleted.
        self.index -= 1
        return modified_strand

    def inject_insertion(self) -> str:
        """
        Inject insertion to the given strand, starting from the base in `index` location of the strand.
        Returns a strand with the injected error.
        :return: a strand with the injected insertion.
        """
        base_insertion_rates = {'A': self.base_error_rates['A']['i'],
                                'T': self.base_error_rates['T']['i'],
                                'C': self.base_error_rates['C']['i'],
                                'G': self.base_error_rates['G']['i']}
        options = list(base_insertion_rates.keys())
        rates = list(base_insertion_rates.values())
        draw = random.choices(options, weights=rates, k=1)
        modified_strand = self.strand[:self.index] + draw[0] + self.strand[self.index:]
        # increment index to approach next original base:
        self.index += 1
        return modified_strand

    def inject_substitution(self) -> str:
        """
        Inject substitution to the given strand, starting from the base in `index` location of the strand.
        Returns a strand with the injected error.
        :return: a strand with the injected substitution.
        """
        base = self.strand[self.index]
        modified_strand = list(self.strand)
        bases = ['A', 'T', 'G', 'C']
        options = []
        for b in bases:
            if b != base:
                options.append(b)
        # Note: 'options' is defined by 'bases' so the order is always the same as in 'bases'.
        # Set rates according to the base:
        rates = [1, 1, 1]
        # if modified_strand[index] == 'G':
        #     rates = []
        draw = random.choices(options, weights=rates, k=1)
        modified_strand[self.index] = draw[0]
        modified_strand = ''.join(modified_strand)
        return modified_strand

    def inject_error(self, error_type: str) -> str:
        """
        Inject the error type to the given strand, starting from the base in `index` location of the strand.
        Returns a strand with the injected error.
        :param error_type: Error type to inject ('d', 'ld', 's', 'pi' as documented)
        :return: a strand with the injected error.
        """
        # check error type and act accordingly:
        if error_type == 'd' or error_type == 'ld':
            return self.inject_deletion(error_type)
        elif error_type == 'pi':  # pre insertion rates are rates for insertion error
            return self.inject_insertion()
        elif error_type == 's':
            return self.inject_substitution()
