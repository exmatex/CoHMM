DistRedis
=========
Distributed Redis implementation is currently on hold pending the authors' ability to figure out a good solution

Problems
=========

As it stands, there are three key problems

1. The older method, implemented with bad, isn't actually functional. A redis slave cannot write and can only read

2.  redis-cluster only distributes the values based upon the keys. This is, in and of itself, a good thing. But it does not redirect gets or sets on its own.
This is problematic as it requires our shmem interface to handle redirects which is inefficient

3. The only solution the authors can think of is to make an intermediary server to do redirects, and we don't want to do that just yet.

Current Plan
=========
Focus on other work and don't deal with this for now. Also consider using the java-based thing with a redis interface.
