Compiling and using .dax reader/writer plugins with ImageJ

To use the plugins:
Copy Dax_Reader.jar and Dax_Writer.jar to your Fiji.app/plugins directory. These can be run from the Plugins menu.


To add support for using File>Open and dragging and dropping, add the modified HandleExtraFileTypes.class to plugins/IO_*.jar
> jar uf /Applications/Fiji.app/plugins/IO_-*.jar HandleExtraFileTypes.class


To compile plugins from source:
> javac -cp /Applications/Fiji.app/jars/ij-*.jar Dax*.java (substituting for the location of your ImageJ jar).
> jar cf Dax_Reader.jar Dax_Reader* HandleExtraFileTypes.java
> jar cf Dax_Writer.jar Dax_Writer* plugins.config 

