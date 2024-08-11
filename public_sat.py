import requests
import json
configUsr = "vishal.b.kolla@gmail.com"
configPwd = "Vishi030!241515"
siteCred = {'identity': configUsr, 'password': configPwd}
uriBase = "https://www.space-track.org"
requestLogin = "/ajaxauth/login"
satellites = "/basicspacedata/query/class/gp/PERIOD/%3C128/orderby/EPOCH%20desc/limit/1000/emptyresult/show/format/json"

with requests.Session() as session:
    resp = session.post(uriBase + requestLogin, data = siteCred)
    if resp.status_code != 200:
        print(resp, "POST fail on login")
    resp = session.get(uriBase + satellites)
    if resp.status_code != 200:
        print(resp)
        session.close()
        exit()
    with open('satellites.json', 'w') as json_file:
        json.dump(resp.json(), json_file, indent=4)
    print("Data saved to satellites.json")
    
    session.close()