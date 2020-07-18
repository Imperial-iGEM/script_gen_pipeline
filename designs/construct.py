# /**
#  * @author Olivia Gallupová
#  * @email olivia.gallupova@gmail.com
#  * @create date 2020-07-17 11:18:14
#  * @modify date 2020-07-17 11:18:14
#  * @desc [description]
#  */


class Construct():
    def __init__(self, input_construct_path, output_sources_paths):

        self.id = ""
        self.parts = ""

        self.input_construct_path = input_construct_path
        self.output_sources_paths = output_sources_paths

        pass

    def get_all_modules(self):
        """ Return all the parts within the final construct
        that are building blocks for assembly """
        raise NotImplementedError
