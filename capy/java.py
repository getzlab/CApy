import jpype

# XXX: this is totally nonfunctional, since jpype doesn't allow for the JVM
#      to be repeatedly started and stopped.
#      can we find a package that actually allows for this?

class javaclass:
	def __init__(self, classname):
		# TODO: add ability to add to this classpath
		jpype.addClassPath("/mnt/j/local/lib/java/*")
		self.classname = classname

	def __enter__(self):
		jpype.startJVM()

		self.inst = jpype.JClass(self.classname);

		return self

	def __exit__(self, *args):
		jpype.shutdownJVM()

def jpa(a):
	return jpype.JArray(jpype.JInt)(a.values)
