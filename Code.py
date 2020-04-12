import os
import datetime
import re

def ParseTweet(tweet):
    return re.findall(r'@(\w+)', tweet)

def ReadTweet(f):
    f.readline()
    TimeLine = f.readline()
    UserLine = f.readline()
    TweetLine = f.readline()

    Timestamp = datetime.datetime.strptime(TimeLine.split('\t')[1].strip(), '%Y-%m-%d %H:%M:%S')
    Username = UserLine.split('\t')[1].strip().split('/')[-1]
    Tweet = TweetLine.split('\t')[1].strip()

    Mentions = ParseTweet(Tweet)

    return (Timestamp, Username, Mentions)

AggregateData = dict()

with open('tweets2009-07.txt', 'r', encoding='utf8') as f:
    i = 0
    while True:
        i += 1
        if i % 100000 == 0:
            print('Processed {} records'.format(i))

        try:
            Timestamp, Username, Mentions = ReadTweet(f)
        except: # End of file
            break
            
        if Timestamp < datetime.datetime(2009, 7, 1, 0, 0, 0) or Timestamp > datetime.datetime(2009, 7, 5, 23, 59, 59):
            continue
        if Timestamp.date() not in AggregateData.keys():
            AggregateData[Timestamp.date()] = dict()

        for mention in Mentions:
            AggregateData[Timestamp.date()][(Username, mention)] = AggregateData[Timestamp.date()].get((Username, mention), 0) + 1

for timestamp, data in AggregateData.items():
    with open(timestamp.strftime('%Y.%m.%d') + '.csv', 'w', encoding = 'utf8') as f:
        f.write('from,to,weight\n')
        for pair, weight in data.items():
            f.write('{},{},{}\n'.format(pair[0], pair[1], weight))
