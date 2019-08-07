





from Communication import Communication
from Test2_worker import worker
from Test2_host import host

    
# put in class and set tags in class as well as comm
# and cell id for host

commun = Communication()   
if commun.get_my_no() == commun.get_host_no():
    host(commun)
else:
    worker(commun)