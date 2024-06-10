def workspace():
    print("location 1: ")
    s = 0
    print("location 2: " + str(s))
    while True:
        print("location 3: " + str(s))
        try:
            print("location 4: " + str(s))
            s= 1/0
            print("location 5: " + str(s))
            break
        except Exception as e:
            print (e)
            print("location 6 " + str(s))
            break
            
    print("location 7: " + str(s))
    
    
if __name__ == '__main__':
    workspace()
