# /**
#  * @author Olivia Gallupová
#  * @email olivia.gallupova@gmail.com
#  * @create date 2020-07-17 11:18:14
#  * @modify date 2020-07-17 11:18:14
#  * @desc [description]
#  */


class Construct():
    def __init__(self, sbol_input):

        self.id = ""
        self.parts = ""

        pass

    def get_all_modules(self):
        """ Return all the parts within the final construct
        that are building blocks for assembly """
        raise NotImplementedError
